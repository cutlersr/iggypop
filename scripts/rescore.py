import pandas as pd
import numpy as np
import argparse
import os
from concurrent.futures import ProcessPoolExecutor

# Purpose of the script:
# This script calculates fidelity scores for a set of overhang sequences
# found in a specified column of an input Excel file. It uses various
# datasets from Potapov et al. and others to compute these scores. 
# The script allows the user to specify additional scoring files, the path
# to the required datasets, and the column name where the overhang sequences
# are located. The results are saved to an output Excel file.

def reverse_complement(seq: str) -> str:
    """
    Generate the reverse complement of a DNA sequence.
    
    Parameters:
    seq (str): The input DNA sequence.

    Returns:
    str: The reverse complement of the input sequence.
    """
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement_dict[base] for base in reversed(seq))

def compute_target_scores(
    seq: str, combined_sequences: list, potapov_data: pd.DataFrame
) -> tuple:
    """
    Compute the on-target and off-target scores for a given sequence based 
    on the provided scoring data.
    
    Parameters:
    seq (str): The sequence to score.
    combined_sequences (list): The list of all sequences and their reverse 
                                complements.
    potapov_data (pd.DataFrame): The dataset containing scoring information.

    Returns:
    tuple: A tuple containing the on-target and off-target scores.
    """
    rc = reverse_complement(seq)
    unique_combined_sequences = list(
        set(combined_sequences + [
            reverse_complement(seq) for seq in combined_sequences
        ])
    )
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
    """
    Calculate the overall fidelity score for a set of sequences based on 
    multiple scoring datasets.
    
    Parameters:
    sequences (list): List of sequences to calculate fidelity scores for.
    potapov_data (pd.DataFrame): The dataset containing scoring information.

    Returns:
    float: The final fidelity score for the set of sequences.
    """
    results = {
        "Sequence": [], "OnTarget": [], "OffTarget": [], 
        "IndividualFidelity": []
    }
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

def process_sequences(sequences, potapov_datasets):
    """
    Process a set of sequences and calculate fidelity scores across 
    multiple datasets.
    
    Parameters:
    sequences (list): List of sequences to process.
    potapov_datasets (list of pd.DataFrame): List of datasets to use for 
                                             scoring.

    Returns:
    tuple: A tuple containing the sequence count and a list of fidelity 
           scores across the datasets.
    """
    sequence_count = len(sequences)
    sequences = [
        seq.strip().upper() for seq in sequences 
        if all(c in 'ATCG' for c in seq.strip().upper())
    ]
    scores = [
        calculate_fidelity_score(sequences, data) 
        for data in potapov_datasets
    ]
    return sequence_count, scores

def parallel_fidelity_scores(sets_of_sequences, potapov_datasets):
    """
    Calculate fidelity scores for multiple sets of sequences in parallel.
    
    Parameters:
    sets_of_sequences (iterable): Iterable of sequence sets to score.
    potapov_datasets (list of pd.DataFrame): List of datasets to use for 
                                             scoring.

    Returns:
    list: A list of tuples containing sequence counts and fidelity scores.
    """
    results = []
    with ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(process_sequences, sequences, potapov_datasets) 
            for sequences in sets_of_sequences
        ]
        for future in futures:
            results.append(future.result())
    return results

def main(input_file, output_file, data_folder, field_name, additional_files):
    """
    Main function to process the input file, calculate fidelity scores, 
    and save the results.
    
    Parameters:
    input_file (str): Path to the input Excel file containing sequence data.
    output_file (str): Path to the output Excel file where results will be 
                       saved.
    data_folder (str): Directory containing the required Potapov datasets.
    field_name (str): Name of the column containing the sets to be scored.
    additional_files (list): List of paths to additional scoring files 
                             provided by the user.
    """
    # Check if the input file exists
    if not os.path.exists(input_file):
        print(
            f"Error: The input file '{input_file}' does not exist. Please "
            "check the file path and try again."
        )
        return

    # Check if the data folder exists
    if not os.path.exists(data_folder):
        print(
            f"Error: The data folder '{data_folder}' does not exist. Please "
            "check the path and try again."
        )
        return

    # Load the input Excel file
    all_data = pd.read_excel(input_file)
    
    # Check if the specified field/column exists in the input file
    if field_name not in all_data.columns:
        print(
            f"Error: The column '{field_name}' does not exist in the input "
            "file."
        )
        if field_name == 'setA':
            print("Tip: You can specify a different column using the '--field' option.")
        return

    # Construct full paths to the Potapov datasets based on the data folder provided
    potapov_datasets = []
    required_files = [
        'FileS04_T4_18h_37C.xlsx',
        'FileS03_T4_18h_25C.xlsx',
        'BsmBI-HFv2_T4.xlsx',
        'BsaI-HFv2_T4.xlsx',
        'FileS06_T7_18h_25C.csv',
        'FileS08_T7_18h_37C.csv'
    ]
    
    for file_name in required_files:
        file_path = os.path.join(data_folder, file_name)
        if not os.path.exists(file_path):
            print(
                f"Error: The required file '{file_name}' was not found in the "
                f"data folder '{data_folder}'."
            )
            return
        if file_name.endswith('.xlsx'):
            potapov_datasets.append(pd.read_excel(file_path))
        elif file_name.endswith('.csv'):
            potapov_datasets.append(pd.read_csv(file_path))
    
    # Load any additional scoring files provided by the user
    for file in additional_files:
        if not os.path.exists(file):
            print(
                f"Warning: The additional scoring file '{file}' does not "
                "exist and will be skipped."
            )
            continue
        if file.endswith('.xlsx'):
            potapov_datasets.append(pd.read_excel(file))
        elif file.endswith('.csv'):
            potapov_datasets.append(pd.read_csv(file))
    
    # Extract sequences from the specified column
    setA_sequences = all_data[field_name].dropna().tolist()
    setA_sets = [sequences.split(',') for sequences in setA_sequences]
    
    # Calculate scores in parallel
    scores = parallel_fidelity_scores(setA_sets, potapov_datasets)
    
    # Define the specific score field names based on the Potapov datasets
    score_fields = [
        'fid_T4_18h_37c', 'fid_T4_18h_25c', 'fid_T4_BsmBI', 'fid_T4_BsaI', 
        'fid_T7_18h_25c', 'fid_T7_18h_37c'
    ]

    # Calculate the size of each set
    all_data['set_size'] = [len(sequences) for sequences in setA_sets]
    all_data['set'] = setA_sequences

    # Append results to the DataFrame with the specific field names
    for index, field in enumerate(score_fields):
        all_data[field] = [score[1][index] for score in scores]

    # Set 'hf_oh_set' as a copy of 'set'
    all_data['hf_oh_set'] = all_data['set']

    # Add 'external_overhangs' if it exists in the input
    if 'external_overhangs' in all_data.columns:
        all_data['external_overhangs'] = all_data['external_overhangs']

    # Ensure the columns are ordered as required
    required_columns = [
        'set_size', 'fidelity', 'set', 'fid_T4_18h_37c', 'fid_T4_18h_25c', 
        'fid_T4_BsmBI', 'fid_T4_BsaI', 'fid_T7_18h_25c', 'fid_T7_18h_37c', 
        'external_overhangs', 'hf_oh_set'
    ]
    final_columns = [
        col for col in required_columns if col in all_data.columns
    ]
    all_data = all_data[final_columns]
    
    # Save the updated DataFrame to the output file
    all_data.to_excel(output_file, index=False)
    
    print(
        f"Excel file has been updated with new fidelity scores and saved to "
        f"{output_file}."
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate fidelity scores for overhang sequences using "
                    "various datasets."
    )
    parser.add_argument(
        '--i', type=str, required=True, 
        help="Input Excel file containing sequence data"
    )
    parser.add_argument(
        '--o', type=str, help="Output Excel file to save the results "
                              "(default: same as input file)"
    )
    parser.add_argument(
        '--d', type=str, default='../../data', 
        help="Directory containing the required Potapov datasets "
             "(default: '../../data')"
    )
    parser.add_argument(
        '--field', type=str, default='setA', 
        help="Name of the column containing the sets to be scored "
             "(default: 'setA')"
    )
    parser.add_argument(
        '--additional_files', nargs='*', 
        help="Additional scoring files (xlsx or csv) to use in scoring"
    )
    
    args = parser.parse_args()
    
    # Determine the output file
    output_file = args.o if args.o else args.i
    
    # Call the main function with the provided arguments
    main(args.i, output_file, args.d, args.field, args.additional_files or [])
