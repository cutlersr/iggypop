import pandas as pd
import numpy as np
from itertools import combinations
from scipy.spatial.distance import pdist, squareform
import argparse

"""
 This script processes an Excel file containing overhang sets and their 
 associated fidelity scores. It filters the data for the top specified 
 percentile of fidelity scores, selects a subset of maximally different 
 sequences based on a Jaccard distance matrix, and outputs the results to 
 an Excel file. The script also generates a distance report that compares 
 the average distances before and after the selection process.

 Expected Input Format:
 - **Excel File**: The file should contain columns for sequence sets and 
   various fidelity scores. The expected columns are:
   - 'set': A string representing the set of sequences, with sequences 
            separated by commas.
   - 'set_size': The number of sequences in the set (this will be calculated 
                 if not present).
   - 'raw_fidelity' or 'fidelity': The raw fidelity score.
   - 'fid_T4_18h_37c', 'fid_T4_18h_25c', 'fid_T4_BsmBI', 'fid_T4_BsaI', 
     'fid_T7_18h_25c', 'fid_T7_18h_37c': Various fidelity scores under 
     different conditions.

 The script outputs:
 - An Excel file containing the selected sets of sequences, sorted by set 
   size and fidelity.
 - A distance report printed to the console that details the average 
   distances for the full and selected subsets.
"""

def read_data(file_path):
    df = pd.read_excel(file_path)
    cols_to_numeric = [
        'set_size', 'fid_T4_18h_37c', 'fid_T4_18h_25c', 
        'fid_T4_BsmBI', 'fid_T4_BsaI', 'fid_T7_18h_25c', 'fid_T7_18h_37c'
    ]
    for col in cols_to_numeric:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    return df

def filter_by_percentile(df, percentile, multi):
    cols_to_filter = [
        'fidelity', 'fid_T4_18h_37c', 'fid_T4_18h_25c', 'fid_T4_BsmBI'
    ]

    def filter_group(group):
        if multi:
            for col in cols_to_filter:
                if col not in group.columns:
                    raise KeyError(f"Column '{col}' not found in the DataFrame.")
                threshold = group[col].quantile(percentile / 100)
                group = group[group[col] >= threshold]
        else:
            if 'fidelity' not in group.columns:
                raise KeyError("Column 'fidelity' not found in the DataFrame.")
            threshold = group['fidelity'].quantile(percentile / 100)
            group = group[group['fidelity'] >= threshold]
        
        top_5 = group.nlargest(5, 'fidelity')
        return pd.concat([group, top_5]).drop_duplicates()

    df_filtered = df.groupby(
        'set_size', group_keys=False
    ).apply(filter_group, include_groups=False)

    df_filtered = df_filtered.drop_duplicates(subset=['set'])
    
    return df_filtered.copy()

def jaccard_distance(set1, set2):
    set1 = set(set1.split(', '))
    set2 = set(set2.split(', '))
    return 1 - len(set1 & set2) / len(set1 | set2)

def create_distance_matrix(sets):
    unique_elements = sorted(
        set.union(*[set(s.split(', ')) for s in sets])
    )
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

def main(file_path, percentile, m, q, multi, output_file):
    df = read_data(file_path)
    df.rename(columns={'raw_fidelity': 'fidelity', 'setA': 'set'}, inplace=True)
    df = df.drop_duplicates(subset=['set'])  

    filtered_df = filter_by_percentile(df, percentile, multi)
    
    if len(filtered_df) < 2 * m:
        print(
            f"Warning: The number of filtered members ({len(filtered_df)}) is "
            f"less than 2*m ({2*m})."
        )
    
    filtered_df['set_size'] = filtered_df['set'].apply(
        lambda x: len(x.split(', '))
    )
    unique_lengths = filtered_df['set_size'].unique()

    selected_rows = pd.DataFrame()
    distance_report = []

    for length in unique_lengths:
        group_df = filtered_df[filtered_df['set_size'] == length]
        top_5_group = group_df.nlargest(5, 'fidelity')
        sets = group_df['set'].tolist()
        distance_matrix = create_distance_matrix(sets)
        
        num_to_select = min(m, len(sets))
        sorted_indices_by_fidelity = np.argsort(
            -group_df['fidelity'].to_numpy()
        )
        initial_indices = sorted_indices_by_fidelity[:q].tolist()

        if len(initial_indices) < num_to_select:
            selected_indices = maxmin_algorithm(
                distance_matrix, num_to_select, initial_indices
            )
        else:
            selected_indices = initial_indices

        selected_group = group_df.iloc[selected_indices]
        selected_group = pd.concat([selected_group, top_5_group]).drop_duplicates()
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

    selected_rows.sort_values(
        by=['set_size', 'fidelity'], ascending=[True, False], inplace=True
    )

    selected_rows.sort_values(
        by=['fid_T4_18h_25c', 'set_size'], ascending=[False, True], inplace=True
    )

    selected_rows.reset_index(drop=True, inplace=True)
    to_remove = []
    for i in range(len(selected_rows) - 1):
        if selected_rows.at[i, 'set_size'] > selected_rows.at[i + 1, 'set_size']:
            to_remove.append(i)
    selected_rows.drop(to_remove, inplace=True)

    columns_to_keep = [
        'set_size', 'fidelity', 'set', 'fid_T4_18h_37c', 'fid_T4_18h_25c', 
        'fid_T4_BsmBI', 'fid_T4_BsaI', 'fid_T7_18h_25c', 'fid_T7_18h_37c'
    ]
    selected_rows = selected_rows[columns_to_keep]

    distance_report_df = pd.DataFrame(distance_report)
    print("\nDistance Report:")
    print(distance_report_df)

    selected_rows.to_excel(output_file, index=False)
    pd.set_option('display.max_rows', None)  
    print("\nSelected Rows:")
    print(selected_rows)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process and thin data based on fidelity and sequence."
    )
    parser.add_argument(
        '--file_path', type=str, required=True, 
        help='Path to the Excel file.'
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
        '--o', type=str, default='hf_ohsets.xlsx', 
        help='Output file name for the final selected sets.'
    )

    args = parser.parse_args()
    main(args.file_path, args.percentile, args.m, args.q, args.multi, args.o)
