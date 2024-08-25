"""
This script processes and filters overhang sequence data from Excel files
to identify maximally different subsets based on fidelity scores. The script
collates data from multiple Excel files in a specified directory, rescoring
the sequences using datasets from Potapov et al. and others. It applies a
percentile cutoff to retain high-fidelity sequences and then selects maximally
distant subsets using the MaxMin algorithm. The final selected sets and a
distance report are saved to an output Excel file.

Usage:
    - Specify the directory containing the input Excel files.
    - Define the percentile cutoff for filtering.
    - Set the number of sequences to select and other parameters as needed.
"""


import os
import pandas as pd
import numpy as np
from itertools import combinations
from scipy.spatial.distance import pdist, squareform
from concurrent.futures import ProcessPoolExecutor
import argparse

def collate_data(directory):
    all_data = []
    for filename in os.listdir(directory):
        if filename.endswith(".xlsx") and not filename.startswith("~$"):
            file_path = os.path.join(directory, filename)
            try:
                param_df = pd.read_excel(file_path, sheet_name='Parameters')
                result_df = pd.read_excel(file_path, sheet_name='Results')
                
                # Rename 'hf_oh_set' to 'set'
                result_df.rename(columns={'hf_oh_set': 'set'}, inplace=True)

                # Add parameter data to each result row
                for param_col in param_df.columns:
                    result_df[param_col] = param_df[param_col].iloc[0]

                # Add set_size column
                result_df['set_size'] = result_df['set'].apply(
                    lambda x: len(x.split(', '))
                )

                all_data.append(result_df)
            except Exception as e:
                print(f"Error processing file {filename}: {e}")
    
    combined_df = pd.concat(all_data, ignore_index=True)
    return combined_df

def read_data(file_path):
    df = pd.read_excel(file_path)
    return df

def filter_by_percentile(df, percentile, multi):
    cols_to_filter = [
        'fidelity', 'fid_T4_18h_37c', 'fid_T4_18h_25c', 'fid_T4_BsmBI'
    ]
    if multi:
        for col in cols_to_filter:
            if col not in df.columns:
                raise KeyError(f"Column '{col}' not found in the DataFrame.")
            df[col] = pd.to_numeric(df[col], errors='coerce')
            threshold = df[col].quantile(percentile / 100)
            df = df[df[col] >= threshold]
    else:
        if 'fidelity' not in df.columns:
            raise KeyError("Column 'fidelity' not found in the DataFrame.")
        df['fidelity'] = pd.to_numeric(df['fidelity'], errors='coerce')
        threshold = df['fidelity'].quantile(percentile / 100)
        df = df[df['fidelity'] >= threshold]

    # Remove identical sets
    df = df.drop_duplicates(subset=['set'])
    return df.copy()

def jaccard_distance(set1, set2):
    set1 = set(set1.split(', '))
    set2 = set(set2.split(', '))
    return 1 - len(set1 & set2) / len(set1 | set2)

def create_distance_matrix(sets):
    unique_elements = sorted(set.union(*[set(s.split(', ')) for s in sets]))
    binary_matrix = np.array([
        [1 if elem in set(s.split(', ')) else 0 for elem in unique_elements] 
        for s in sets
    ])
    distance_matrix = pdist(binary_matrix, metric='jaccard')
    return squareform(distance_matrix)

def maxmin_algorithm(distance_matrix, m, initial_indices):
    selected_indices = initial_indices.copy()
    remaining_indices = set(range(len(distance_matrix))) - set(selected_indices)

    while len(selected_indices) < m and remaining_indices:
        max_distance = -1
        next_index = -1
        for index in remaining_indices:
            min_distance_to_selected = min(
                distance_matrix[index, selected] 
                for selected in selected_indices
            )
            if min_distance_to_selected > max_distance:
                max_distance = min_distance_to_selected
                next_index = index
        remaining_indices.remove(next_index)
        selected_indices.append(next_index)
    
    return selected_indices

def reverse_complement(seq: str) -> str:
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement_dict[base] for base in reversed(seq))

def compute_target_scores(
    seq: str, combined_sequences: list, potapov_data: pd.DataFrame
) -> tuple:
    rc = reverse_complement(seq)
    unique_combined_sequences = list(set(
        combined_sequences + [reverse_complement(seq) for seq in combined_sequences]
    ))
    rc_rows = potapov_data[potapov_data["Overhang"] == rc]
    seq_rows = potapov_data[potapov_data["Overhang"] == seq]
    on_target_score = rc_rows[seq].values[0] + seq_rows[rc].values[0]
    off_target_score = rc_rows[unique_combined_sequences].drop(
        columns=[seq]
    ).sum().sum() + seq_rows[unique_combined_sequences].drop(
        columns=[rc]
    ).sum().sum()
    return on_target_score, off_target_score

def calculate_fidelity_score(
    sequences: list, potapov_data: pd.DataFrame
) -> float:
    results = {
        "Sequence": [], "OnTarget": [], "OffTarget": [], 
        "IndividualFidelity": []
    }
    sequences = list(sequences)
    rc_sequences = [reverse_complement(seq) for seq in sequences]
    combined_sequences = list(set(sequences + rc_sequences))
    for seq in sequences:
        on_target_score, off_target_score = compute_target_scores(
            seq, combined_sequences, potapov_data
        )
        individual_fidelity = 1 - (
            off_target_score / (on_target_score + off_target_score)
        )
        results["Sequence"].append(seq)
        results["OnTarget"].append(on_target_score)
        results["OffTarget"].append(off_target_score)
        results["IndividualFidelity"].append(individual_fidelity)
    final_fidelity = np.prod(results["IndividualFidelity"])
    return final_fidelity

def process_sequences(sequences, potapov_data_list):
    sequence_count = len(sequences)
    sequences = [
        seq.strip().upper() for seq in sequences 
        if all(c in 'ATCG' for c in seq.strip().upper())
    ]
    fidelity_scores = [
        calculate_fidelity_score(sequences, potapov_data) 
        for potapov_data in potapov_data_list
    ]
    return sequence_count, fidelity_scores

def parallel_fidelity_scores(sets_of_sequences, potapov_data_list):
    results = []
    with ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(process_sequences, sequences, potapov_data_list) 
            for sequences in sets_of_sequences
        ]
        for future in futures:
            sequence_count, scores = future.result()
            results.append((sequence_count, scores))
    return results

def main(directory, percentile, m, q, multi, output_file):
    combined_df = collate_data(directory)
    print("Combined DataFrame after collating:")
    print(combined_df.head())

    # Write combined data to a temporary file
    temp_file = "temp_collated_data.xlsx"
    combined_df.to_excel(temp_file, index=False)
    
    # Load the Potapov data files
    potapov_data_files = [
        '../data/FileS04_T4_18h_37C.xlsx',
        '../data/FileS03_T4_18h_25C.xlsx',
        '../data/BsmBI-HFv2_T4.xlsx',
        '../data/BsaI-HFv2_T4.xlsx',
        '../data/FileS06_T7_18h_25C.csv',
        '../data/FileS08_T7_18h_37C.csv'
    ]
    potapov_data_list = [
        pd.read_excel(f) if f.endswith('.xlsx') else pd.read_csv(f) 
        for f in potapov_data_files
    ]

    # Rescore the sequences
    combined_df = pd.read_excel(temp_file)
    setA_sequences = combined_df['set'].dropna().unique().tolist()
    setA_sets = combined_df['set'].dropna().apply(lambda x: x.split(','))

    # Calculate scores in parallel
    scores = parallel_fidelity_scores(setA_sets, potapov_data_list)
    
    # Ensure the scores are extracted from the tuples correctly
    score_fields = [
        'fid_T4_18h_37c', 'fid_T4_18h_25c', 'fid_T4_BsmBI', 'fid_T4_BsaI', 
        'fid_T7_18h_25c', 'fid_T7_18h_37c'
    ]
    for index, field in enumerate(score_fields):
        combined_df[field] = [score[1][index] for score in scores]

    # Add set_size column
    combined_df['set_size'] = combined_df['set'].apply(
        lambda x: len(x.split(', '))
    )

    # Save the rescored data
    rescored_file = "rescored_data.xlsx"
    combined_df.to_excel(rescored_file, index=False)
    print("Combined DataFrame after rescoring:")
    print(combined_df.head())
    
    # Load rescored data
    rescored_df = pd.read_excel(rescored_file)

    print(rescored_df)

    # Filter by percentile
    filtered_df = filter_by_percentile(rescored_df, percentile, multi)
    print("Filtered DataFrame:")
    print(filtered_df.head())
    
    if len(filtered_df) < 2 * m:
        print(
            f"Warning: The number of filtered members ({len(filtered_df)}) is "
            f"less than 2*m ({2*m})."
        )
    
    # Group by set length
    unique_lengths = filtered_df['set_size'].unique()

    print(unique_lengths)

    selected_rows = pd.DataFrame()

    print(selected_rows)

    distance_report = []

    for length in unique_lengths:
        group_df = filtered_df[filtered_df['set_size'] == length]
        sets = group_df['set'].tolist()
        distance_matrix = create_distance_matrix(sets)

        print(distance_matrix)
        
        # Handle cases where the number of available sets is less than m
        num_to_select = min(m, len(sets))
        
        # Include the top q highest-scoring members
        sorted_indices_by_fidelity = np.argsort(-group_df['fidelity'].to_numpy())
        initial_indices = sorted_indices_by_fidelity[:q].tolist()

        if len(initial_indices) < num_to_select:
            # Apply maxmin algorithm starting from the initial indices
            selected_indices = maxmin_algorithm(
                distance_matrix, num_to_select, initial_indices
            )
        else:
            selected_indices = initial_indices

        selected_group = group_df.iloc[selected_indices]
        selected_rows = pd.concat([selected_rows, selected_group])

        avg_distance_filtered = np.mean(distance_matrix)
        avg_distance_selected = np.mean(
            distance_matrix[selected_indices][:, selected_indices].flatten()
        )

        distance_report.append({
            'set_size': length,
            'filtered': avg_distance_filtered,
            'final': avg_distance_selected
        })

    # Sort the final DataFrame by 'set_size' ascending and 'fid_T4_18h_25c' descending
    selected_rows.sort_values(
        by=['set_size', 'fid_T4_18h_25c'], ascending=[True, False], inplace=True
    )
    print("Selected Rows DataFrame after sorting:")
    print(selected_rows.head())

    # Select only the desired columns
    columns_to_keep = [
        'set_size', 'fidelity', 'set', 'fid_T4_18h_37c', 'fid_T4_18h_25c', 
        'fid_T4_BsmBI', 'fid_T4_BsaI', 'fid_T7_18h_25c', 'fid_T7_18h_37c'
    ]
    selected_rows = selected_rows[columns_to_keep]

    # Create distance report DataFrame
    distance_report_df = pd.DataFrame(distance_report)

    # Print distance report
    print("\nDistance Report:")
    print(distance_report_df)

    # Write the final selected sets to an Excel file
    selected_rows.to_excel(output_file, index=False)

    pd.set_option('display.max_rows', None)  # Ensure the full DataFrame is printed
    print("\nSelected Rows:")
    print(selected_rows)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process and thin data based on fidelity and sequence."
    )
    parser.add_argument(
        '--directory', type=str, required=True, 
        help='Directory containing the Excel files.'
    )
    parser.add_argument(
        '--percentile', type=float, required=True, 
        help='Percentile cutoff for filtering fidelity.'
    )
    parser.add_argument(
        '--m', type=int, required=True, 
        help='Number of maximally different members to select for each set length.'
    )
    parser.add_argument(
        '--q', type=int, default=2, 
        help='Number of highest-scoring members to include initially.'
    )
    parser.add_argument(
        '--multi', action='store_true', 
        help='Apply the percentile cutoff to multiple columns.'
    )
    parser.add_argument(
        '--output_file', type=str, default='hf_ohsets.xlsx', 
        help='Output file name for the final selected sets.'
    )

    args = parser.parse_args()
    main(
        args.directory, args.percentile, args.m, args.q, 
        args.multi, args.output_file
    )
