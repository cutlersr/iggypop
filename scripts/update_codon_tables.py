import os
import pandas as pd

def calculate_relative_frequencies(row):
    total_codons = sum([row[codon] for codon in codon_list])
    return {codon: row[codon] / total_codons for codon in codon_list}

def main():
    # Define the target directory
    target_directory = '.venv/lib/python3.12/site-packages/codon_usage_data/tables'
    os.makedirs(target_directory, exist_ok=True)  # Ensure the directory exists
    
    # Read the master data file
    master_data = pd.read_csv('path_to_master_file.csv', sep='\t')
    
    # List of codons for each amino acid
    codon_list = [
        # Add your complete list of codons here
    ]
    
    amino_acid_to_codons = {
        # Define your mapping of amino acids to codons
    }
    
    # Generate and save individual codon usage files
    for index, row in master_data.iterrows():
        species_name = row['Species'].replace(' ', '_')
        taxid = row['Taxid']
        filename = f"{species_name}_{taxid}.csv"
        file_path = os.path.join(target_directory, filename)  # Full path to the file
        
        # Create the codon usage table
        frequencies = calculate_relative_frequencies(row)
        data = []
        for aa, codons in amino_acid_to_codons.items():
            for codon in codons:
                data.append([aa, codon, frequencies.get(codon, 0)])
                
        df = pd.DataFrame(data, columns=['amino_acid', 'codon', 'relative_frequency'])
        df.to_csv(file_path, index=False)
        print(f"Saved codon usage table to {file_path}")

if __name__ == "__main__":
    main()

