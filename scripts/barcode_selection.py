import numpy as np
import pandas as pd
import argparse
import re
import editdistance
import os
from datetime import datetime

"""
 This script analyzes a set of candidate barcodes to identify maximally 
 different subsets. It does this by computing pairwise edit distances for 
 all candidates, then selecting maximally distant subsets. It pre-filters 
 sequences with homonucleotide repeats. The script can also save and load 
 distance matrices for efficient reuse if needed.

 The input CSV file needs either or both of the following columns:
 - 'F_seq': Forward sequences of the barcodes.
 - 'R_seq_rc': Reverse complement sequences of different barcodes.

 If both columns are present, the candidates are treated as distinct 
 sequences.

 We used this script to select our barcodes for nanopore amplicon sequencing,
 identifying a subset of 18 18-mers with an average edit distance of ~10
 and minimum of 8.
"""

# Function to compute pairwise edit distances
def compute_edit_distances(sequences):
    n = len(sequences)
    distance_matrix = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            distance = editdistance.eval(sequences[i], sequences[j])
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance
    return distance_matrix

# Function to select a maximally distant subset of sequences
def select_maximally_distant_subset(sequences, distance_matrix, n):
    num_sequences = len(sequences)
    selected_members = np.zeros(n, dtype=int)
    selected_members[0] = np.argmax(np.sum(distance_matrix, axis=1))
    for i in range(1, n):
        min_distances = np.min(
            distance_matrix[:, selected_members[:i]], axis=1
        )
        selected_members[i] = np.argmax(min_distances)
    selected_sequences = [sequences[i] for i in selected_members]
    return selected_sequences, selected_members

# Function to filter out sequences with any single base that repeats three 
# or more times consecutively
def filter_sequences(sequences):
    original_count = len(sequences)
    filtered_sequences = sequences[
        ~sequences.str.contains(r'(A{3,}|T{3,}|C{3,}|G{3,})')
    ]
    filtered_count = len(filtered_sequences)
    print(
        f"Filtered out {original_count - filtered_count} sequences due to "
        "homonucleotide repeats."
    )
    return filtered_sequences

# Function to compute the reverse complement of a DNA sequence
def reverse_complement(seq):
    complement = str.maketrans('ATCG', 'TAGC')
    return seq.translate(complement)[::-1]

# Function to save the distance matrix and sequences to a compressed file
def save_distance_matrix(distance_matrix, sequences, file_name):
    np.savez_compressed(
        file_name, distance_matrix=distance_matrix, sequences=sequences
    )

# Function to load the distance matrix and sequences from a compressed file
def load_distance_matrix_and_sequences(file_name):
    with np.load(file_name) as data:
        if 'distance_matrix' not in data or 'sequences' not in data:
            raise KeyError(
                f"The file {file_name} does not contain both "
                "'distance_matrix' and 'sequences'."
            )
        return data['distance_matrix'], data['sequences']

# Function to print statistics for the selected subset
def print_selected_set_stats(distance_matrix, selected_indices):
    n = len(selected_indices)
    subset_matrix = np.zeros((n, n), dtype=int)
    
    for i in range(n):
        for j in range(n):
            subset_matrix[i, j] = distance_matrix[
                selected_indices[i], selected_indices[j]
            ]
    
    # Calculate statistics on the subset matrix
    distances = subset_matrix[np.triu_indices(n, k=1)]
    
    print("Statistics for the selected subset:")
    print(f"Minimum edit distance: {np.min(distances)}")
    print(f"Maximum edit distance: {np.max(distances)}")
    print(f"Mean edit distance: {np.mean(distances):.2f}")
    print(f"Median edit distance: {np.median(distances)}")

# Setting up command line argument parsing
parser = argparse.ArgumentParser(
    description="Select a maximally distant subset of barcodes."
)
parser.add_argument(
    "--file_path", type=str, help="Path to the file containing the barcodes."
)
parser.add_argument(
    "--num_barcodes", type=int, help="Number of barcodes to select."
)
parser.add_argument(
    "--distance_matrix_file", type=str, default=None, 
    help="Path to a precomputed distance matrix file (.npz)."
)
parser.add_argument(
    "--save_matrix", type=str, default="distance_matrix.npz", 
    help="Path to save the computed distance matrix (.npz)."
)
parser.add_argument(
    "--sample", type=int, default=None, 
    help="Number of random samples to take from the input data for testing."
)
parser.add_argument(
    "--ignore_reverse_complements", action="store_true", 
    help="Ignore reverse complements in the analysis."
)
args = parser.parse_args()

# Track stop time
start_time = datetime.now()
print(f"Script started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

# Check if a precomputed distance matrix is provided
if args.distance_matrix_file and os.path.exists(args.distance_matrix_file):
    print(
        f"Loading distance matrix and sequences from {args.distance_matrix_file}..."
    )
    try:
        distance_matrix, sequences = load_distance_matrix_and_sequences(
            args.distance_matrix_file
        )
        sequences = sequences.tolist()  # Ensure the sequences are in list format
    except KeyError as e:
        print(e)
        print(
            "The distance matrix does not contain the sequences. Please "
            "regenerate the matrix."
        )
        exit(1)
else:
    # Read from CSV file
    df = pd.read_csv(args.file_path, dtype=str)

    # Sample the data if the sample parameter is provided
    if args.sample:
        df = df.sample(n=args.sample, random_state=42).reset_index(drop=True)

    # Convert sequences to strings to handle potential mixed types
    df['F_seq'] = df['F_seq'].astype(str)
    df['R_seq_rc'] = df['R_seq_rc'].astype(str)

    # Filter sequences to remove any with triple or more consecutive 
    # identical letters before generating reverse complements
    df_filtered = df[
        (~df['F_seq'].str.contains(r'(A{3,}|T{3,}|C{3,}|G{3,})')) & 
        (~df['R_seq_rc'].str.contains(r'(A{3,}|T{3,}|C{3,}|G{3,})'))
    ]

    if args.ignore_reverse_complements:
        sequences = pd.concat(
            [df_filtered['F_seq'], df_filtered['R_seq_rc']]
        ).reset_index(drop=True).tolist()
    else:
        sequences = pd.concat([
            df_filtered['F_seq'],
            df_filtered['R_seq_rc'],
            df_filtered['F_seq'].apply(reverse_complement),
            df_filtered['R_seq_rc'].apply(reverse_complement)
        ]).reset_index(drop=True).tolist()

    # Remove any NaN values or invalid sequences before computing the distance matrix
    sequences = [
        seq for seq in sequences if pd.notna(seq) and seq != 'nan'
    ]

    print("Computing the pairwise edit distances...")
    distance_matrix = compute_edit_distances(sequences)
    
    # Save the distance matrix and sequences
    save_distance_matrix(distance_matrix, sequences, args.save_matrix)
    print(f"Distance matrix saved to {args.save_matrix}.")

# Check if the number of barcodes requested is less than the total available
if args.num_barcodes > len(sequences):
    raise ValueError(
        "Requested number of barcodes is greater than the number available."
    )

# Select a subset of maximally distant sequences
selected_sequences, selected_indices = select_maximally_distant_subset(
    sequences, distance_matrix, args.num_barcodes
)

# Print and write the selected sequences
print("Selected sequences:")
with open("barcodes.txt", "w") as file:
    for seq in selected_sequences:
        print(seq)
        file.write(seq + '\n')

# Print statistics for the selected subset
print_selected_set_stats(distance_matrix, selected_indices)

# Track stop time
stop_time = datetime.now()
print(f"Script finished at {stop_time.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Total runtime: {stop_time - start_time}")
