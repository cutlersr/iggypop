import pandas as pd
import random
from scipy.spatial.distance import pdist, squareform
import numpy as np
import argparse

# Set up command-line arguments
parser = argparse.ArgumentParser(description="Generate random mutant sequences and select a diverse subset.")
parser.add_argument('--i', type=str, default="seq_data.tsv")
parser.add_argument('--name', type=str, default="nnLuz")
parser.add_argument('--n', type=int, default=1000, help="Number of random sequences to generate")
parser.add_argument('--min', type=int, default=14, help="Minimum edit distance from the wild type required")
parser.add_argument('--max', type=int, default=20, help="Maximum edit distance from the wild type allowed")
parser.add_argument('--n_seqs', type=int, default=5, help="Number of sequences to return at the end")

args = parser.parse_args()

# Read the data from the file
file_path = args.i  # Replace with your actual file path
data = pd.read_csv(file_path, sep='\t')

# Generate random amino acid sequences
random_sequences = []
for _ in range(args.n):
    sequence = "".join(
        random.choice([wt_aa, alt_aa]) if pd.notna(alt_aa) and alt_aa != "" else wt_aa
        for wt_aa, alt_aa in zip(data['wt_aa'], data['alt_aa'])
    )
    random_sequences.append(sequence)

# Compute pairwise edit distances
wild_type_sequence = "".join(data['wt_aa'])
edit_distances = []
for seq in random_sequences:
    distance = sum(1 for a, b in zip(wild_type_sequence, seq) if a != b)
    edit_distances.append(distance)

# Find sequences with at least the specified minimum distance and at most the specified maximum distance from the wild type
sequences_with_min_max_edit = [seq for seq, dist in zip(random_sequences, edit_distances) if args.min <= dist <= args.max]

# Calculate the pairwise distance matrix for the subset
pairwise_dist_matrix = squareform(pdist(np.array(sequences_with_min_max_edit).reshape(-1, 1),
                                        metric=lambda u, v: sum(1 for x, y in zip(u[0], v[0]) if x != y)))

# Select a maximally diverse subset of sequences using Farthest Point Sampling (FPS)
def farthest_point_sampling(dist_matrix, n_points):
    selected_indices = [0]  # Start with the first sequence
    while len(selected_indices) < n_points and len(selected_indices) < len(dist_matrix):
        max_min_distance = -1
        best_index = -1
        for i in range(len(dist_matrix)):
            if i not in selected_indices:
                min_distance_to_selected = min(dist_matrix[i, j] for j in selected_indices)
                if min_distance_to_selected > max_min_distance:
                    max_min_distance = min_distance_to_selected
                    best_index = i
        if best_index != -1:
            selected_indices.append(best_index)
    return selected_indices

# Select indices using FPS
selected_indices = farthest_point_sampling(pairwise_dist_matrix, args.n_seqs)
selected_sequences = [sequences_with_min_max_edit[i] for i in selected_indices]

# Generate nucleic acid sequences for the selected amino acid sequences
nucleic_acid_sequences = []
for seq in selected_sequences:
    nucleic_seq = "".join(
        data.loc[idx, 'alt_seq'] if aa == data.loc[idx, 'alt_aa'] and pd.notna(data.loc[idx, 'alt_seq']) else data.loc[idx, 'wt_seq']
        for idx, aa in enumerate(seq)
    )
    nucleic_acid_sequences.append(nucleic_seq)

# Write the selected protein sequences to a FASTA file
protein_output_path = 'selected_protein_sequences.fasta'
with open(protein_output_path, 'w') as f:
    for i, aa_seq in enumerate(selected_sequences, start=1):
        f.write(f">{args.name}_protein_v{i}\n{aa_seq}\n")

# Write the selected nucleic acid sequences to a FASTA file
nucleic_output_path = 'selected_nucleic_acid_sequences.fasta'
with open(nucleic_output_path, 'w') as f:
    for i, nuc_seq in enumerate(nucleic_acid_sequences, start=1):
        f.write(f">{args.name}_nucleic_v{i}\n{nuc_seq}\n")

# Write the edit distance matrix for the selected proteins
edit_distance_matrix_path = 'edit_distance_matrix.csv'
np.savetxt(edit_distance_matrix_path, pairwise_dist_matrix[np.ix_(selected_indices, selected_indices)], delimiter=",", fmt="%d")

# Write a table of edit distances to the wild type sequence
edit_distances_to_wt_path = 'edit_distances_to_wild_type.csv'
pd.DataFrame({
    'Sequence': [f'{args.name}_protein_v{i+1}' for i in range(len(selected_sequences))],
    'Edit Distance to WT': [sum(1 for a, b in zip(wild_type_sequence, seq) if a != b) for seq in selected_sequences]
}).to_csv(edit_distances_to_wt_path, index=False)

print(f"Protein sequences written to {protein_output_path}")
print(f"Nucleic acid sequences written to {nucleic_output_path}")
print(f"Edit distance matrix written to {edit_distance_matrix_path}")
print(f"Edit distances to wild type written to {edit_distances_to_wt_path}")
