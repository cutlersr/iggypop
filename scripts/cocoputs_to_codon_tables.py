import os
import pandas as pd

# List of all codons with U instead of T
codon_list = [
    'UUU', 'UUC', 'UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG',
    'AUU', 'AUC', 'AUA', 'AUG', 'GUU', 'GUC', 'GUA', 'GUG',
    'UAU', 'UAC', 'UAA', 'UAG', 'CAU', 'CAC', 'CAA', 'CAG',
    'AAU', 'AAC', 'AAA', 'AAG', 'GAU', 'GAC', 'GAA', 'GAG',
    'UCU', 'UCC', 'UCA', 'UCG', 'CCU', 'CCC', 'CCA', 'CCG',
    'ACU', 'ACC', 'ACA', 'ACG', 'GCU', 'GCC', 'GCA', 'GCG',
    'UGU', 'UGC', 'UGA', 'UGG', 'CGU', 'CGC', 'CGA', 'CGG',
    'AGU', 'AGC', 'AGA', 'AGG', 'GGU', 'GGC', 'GGA', 'GGG'
]

# Mapping from amino acids to codons with U instead of T
amino_acid_to_codons = {
    '*': ['UAA', 'UAG', 'UGA'],  # Stop codons
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'C': ['UGU', 'UGC'],
    'D': ['GAU', 'GAC'],
    'E': ['GAA', 'GAG'],
    'F': ['UUU', 'UUC'],
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],
    'H': ['CAU', 'CAC'],
    'I': ['AUU', 'AUC', 'AUA'],
    'K': ['AAA', 'AAG'],
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'M': ['AUG'],
    'N': ['AAU', 'AAC'],
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'Q': ['CAA', 'CAG'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    'W': ['UGG'],
    'Y': ['UAU', 'UAC']
}

def calculate_relative_frequencies(row):
    total_codons = sum(row[codon] for codon in codon_list if codon in row)
    return {codon: (row[codon] / total_codons if total_codons > 0 else 0) for codon in codon_list if codon in row}

def format_short_name(short_name):
    parts = short_name.split('_')
    return f"{parts[0][0]}_{parts[1]}" if len(parts) > 1 else short_name

def main():
    target_directory = '/usr/local/lib/python3.9/site-packages/codon_usage_data/tables'
    os.makedirs(target_directory, exist_ok=True)
    
    codon_usage_directory = '/usr/local/lib/python3.9/site-packages/codon_usage_data'
    
    master_data = pd.read_csv('data/cleaned_coco.tsv', sep='\t')
    master_data['Total Codons'] = master_data[codon_list].sum(axis=1)
    best_examples = master_data.loc[master_data.groupby('Species')['Total Codons'].idxmax()]

    # Save the thinned master data
    best_examples.to_csv('data/cleaned_coco_thinned.tsv', sep='\t', index=False)
    
    # Save organism data
    organisms_data = best_examples[['short_name', 'Taxid']].copy()
    organisms_data['organism'] = organisms_data['short_name'].apply(format_short_name)
    organisms_data.drop('short_name', axis=1, inplace=True)
    organisms_csv_path = os.path.join(codon_usage_directory, 'organisms.csv')
    organisms_data.to_csv(organisms_csv_path, index=False)
    
    table_count = 0
    for index, row in best_examples.iterrows():
        short_name = format_short_name(row['short_name'])
        taxid = row['Taxid']
        filename = f"{short_name}_{taxid}.csv"
        file_path = os.path.join(target_directory, filename)
        
        frequencies = calculate_relative_frequencies(row)
        data = []
        for aa, codons in amino_acid_to_codons.items():
            for codon in codons:
                data.append([aa, codon, frequencies.get(codon, 0)])
                
        df = pd.DataFrame(data, columns=['amino_acid', 'codon', 'relative_frequency'])
        df.to_csv(file_path, index=False)
        table_count += 1

    print(f"Saved {table_count} codon usage tables to {target_directory}")
    print(f"Organism data saved to {organisms_csv_path}")

if __name__ == "__main__":
    main()
