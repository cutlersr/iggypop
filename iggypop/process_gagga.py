import pandas as pd
import numpy as np
from itertools import combinations
from scipy.spatial.distance import pdist, squareform
import argparse

def read_data(file_path):
    return pd.read_excel(file_path)

def filter_by_percentile(df, percentile):
    threshold = df['raw_fidelity'].quantile(percentile / 100)
    return df[df['raw_fidelity'] >= threshold]

def jaccard_distance(set1, set2):
    set1 = set(set1.split(', '))
    set2 = set(set2.split(', '))
    return 1 - len(set1 & set2) / len(set1 | set2)

def create_distance_matrix(sets):
    distance_matrix = pdist(sets, lambda u, v: jaccard_distance(u, v))
    return squareform(distance_matrix)

def select_maximally_different(distance_matrix, m):
    selected_indices = []
    remaining_indices = set(range(len(distance_matrix)))

    while len(selected_indices) < m and remaining_indices:
        if not selected_indices:
            next_index = remaining_indices.pop()
        else:
            max_distance = -1
            next_index = -1
            for index in remaining_indices:
                min_distance = min(distance_matrix[index, selected] for selected in selected_indices)
                if min_distance > max_distance:
                    max_distance = min_distance
                    next_index = index
            remaining_indices.remove(next_index)
        
        selected_indices.append(next_index)
    
    return selected_indices

def main(file_path, percentile, m):
    df = read_data(file_path)
    filtered_df = filter_by_percentile(df, percentile)
    sets = filtered_df['setA'].tolist()
    distance_matrix = create_distance_matrix(sets)
    selected_indices = select_maximally_different(distance_matrix, m)
    selected_rows = filtered_df.iloc[selected_indices]
    print(selected_rows)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process and thin data based on raw_fidelity and setA.")
    parser.add_argument('file_path', type=str, help='Path to the Excel file.')
    parser.add_argument('percentile', type=float, help='Percentile cutoff for filtering raw_fidelity.')
    parser.add_argument('m', type=int, help='Number of maximally different members to select.')

    args = parser.parse_args()
    main(args.file_path, args.percentile, args.m)
