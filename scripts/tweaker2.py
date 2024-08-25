
"""
This script extracts DNA sequences from an iggypop log file, groups them by 
their accessions, computes pairwise edit distances of the sequences for a 
given accession, and then selects 'n' maximally different subset of sequences 
based on these distances using the maxmin algorithm. The selected sequences 
and the distance matrices are saved to output files.

Expected Input Format:
- **Log File**: iggypop log file with FASTA formatted chiseled sequences

The script outputs:
- A FASTA file containing the selected subset of maximally different sequences.
- A text file containing the full distance matrix and the subset distance matrix 
  for each base gene group, along with the average distances for both the full 
  set and the subset.
"""

import argparse
import re
from itertools import combinations
import numpy as np
import pandas as pd
from Levenshtein import distance as levenshtein_distance
import os

def extract_sequences(file_path, n):
    sequences = []
    with open(file_path, 'r') as file:
        current_seq_id = None
        current_sequence = []
        capture_sequence = False

        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_seq_id:
                    sequences.append((current_seq_id, ''.join(current_sequence)))
                    current_seq_id = None
                    current_sequence = []

                current_seq_id_match = re.match(r'>(.*\.\d+_chisel)', line)
                if current_seq_id_match:
                    current_seq_id = current_seq_id_match.group(1)
                    capture_sequence = True
                else:
                    capture_sequence = False
            elif capture_sequence and re.match(r'^[ATCG]+$', line):
                current_sequence.append(line)

        if current_seq_id:
            sequences.append((current_seq_id, ''.join(current_sequence)))

    return sequences

def compute_distance_matrix(sequences):
    num_sequences = len(sequences)
    matrix = np.zeros((num_sequences, num_sequences))

    for i, j in combinations(range(num_sequences), 2):
        dist = levenshtein_distance(sequences[i][1], sequences[j][1])
        max_len = max(len(sequences[i][1]), len(sequences[j][1]))
        percentage = (dist / max_len) * 100
        matrix[i, j] = percentage
        matrix[j, i] = percentage

    return matrix

def group_sequences_by_base_gene(sequences):
    base_gene_groups = {}
    for seq_id, sequence in sequences:
        base_gene = re.match(r'(.*)\.\d+', seq_id).group(1)
        if base_gene not in base_gene_groups:
            base_gene_groups[base_gene] = []
        base_gene_groups[base_gene].append((seq_id, sequence))

    return base_gene_groups

def maxmin_selection(distance_matrix, n):
    """
    Perform MaxMin selection on a distance matrix to select a diverse subset 
    of elements.

    The MaxMin algorithm is designed to iteratively select elements that are 
    maximally distant from each other, ensuring that the selected subset 
    captures the diversity of the entire dataset. The algorithm selects the 
    pair of elements with the maximum distance between them; subsequent 
    selections are made by choosing the next element that has the maximum 
    minimum distance to any element already in the selected subset.

    Reference:
    Butina, D. (1999). Unsupervised Data Base Clustering Based on Daylight's 
    Fingerprint and Tanimoto Similarity: A Fast and Automated Way To Cluster 
    Small and Large Data Sets. Journal of Chemical Information and Computer 
    Sciences, 39(4), 747-750. doi:10.1021/ci9803381

    Parameters:
    - distance_matrix (ndarray): A square matrix containing the pairwise 
      distances between elements.
    - n (int): The number of elements to select.

    Returns:
    - list: The indices of the selected elements.
    """
    selected_indices = []
    num_sequences = distance_matrix.shape[0]

    # Select the initial pair with the maximum distance
    max_dist = -1
    initial_pair = (0, 1)
    for i, j in combinations(range(num_sequences), 2):
        if distance_matrix[i, j] > max_dist:
            max_dist = distance_matrix[i, j]
            initial_pair = (i, j)

    selected_indices.extend(initial_pair)

    # Select the rest of the members iteratively
    while len(selected_indices) < n:
        max_min_dist = -1
        next_index = -1
        for i in range(num_sequences):
            if i not in selected_indices:
                min_dist = min(
                    distance_matrix[i, selected_index]
                    for selected_index in selected_indices
                )
                if min_dist > max_min_dist:
                    max_min_dist = min_dist
                    next_index = i
        selected_indices.append(next_index)

    return selected_indices

def main():
    parser = argparse.ArgumentParser(
        description='Extract sequences from log file and compute pairwise edit distances.'
    )
    parser.add_argument('--file', required=True, help='Path to the log file')
    parser.add_argument(
        '--tweak_n', default=5, type=int, required=True, help='Number of output sequences'
    )
    
    args = parser.parse_args()
    sequences = extract_sequences(args.file, args.tweak_n)
    
    base_gene_groups = group_sequences_by_base_gene(sequences)
    output_dir = os.path.dirname(args.file)

    with open(os.path.join(output_dir, 'max_diff_subsets.fasta'), 'w') as fasta_file, \
         open(os.path.join(output_dir, 'distances.txt'), 'w') as distance_file:

        for base_gene, seqs in base_gene_groups.items():
            print(f'Base Gene: {base_gene}')
            matrix = compute_distance_matrix(seqs)
            selected_indices = maxmin_selection(
                matrix, min(len(seqs), args.tweak_n)
            )
            selected_sequences = [seqs[i] for i in selected_indices]
            
            # Write selected sequences to FASTA file
            for seq_id, sequence in selected_sequences:
                fasta_file.write(f'>{seq_id}\n{sequence}\n')
            
            # Write full distance matrix to distances.txt
            full_df = pd.DataFrame(matrix, index=[s[0] for s in seqs])
            full_df = full_df.astype(int).astype(str)
            distance_file.write(f'Full Distance Matrix for {base_gene}:\n')
            distance_file.write(full_df.to_string(index=True, header=False))
            distance_file.write('\n\n')
            
            # Write subset distance matrix to distances.txt
            subset_matrix = matrix[np.ix_(selected_indices, selected_indices)]
            subset_df = pd.DataFrame(
                subset_matrix, index=[s[0] for s in selected_sequences]
            )
            subset_df = subset_df.astype(int).astype(str)
            distance_file.write(f'Subset Distance Matrix for {base_gene}:\n')
            distance_file.write(subset_df.to_string(index=True, header=False))
            distance_file.write('\n\n')
            
            # Calculate and write average distance for the entire set
            avg_distance_full = np.mean(matrix[np.triu_indices(len(matrix), 1)])
            distance_file.write(
                f'Average Distance for Full Set of {base_gene}: {avg_distance_full:.2f}%\n'
            )
            
            # Calculate and write average distance for the selected subset
            avg_distance_subset = np.mean(
                subset_matrix[np.triu_indices(len(subset_matrix), 1)]
            )
            distance_file.write(
                f'Average Distance for Selected Subset of {base_gene}: {avg_distance_subset:.2f}%\n'
            )
            distance_file.write('\n\n')

if __name__ == '__main__':
    main()
