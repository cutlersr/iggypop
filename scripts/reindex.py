import pandas as pd
import argparse

def reformat_fragments(xlsx_file, csv_file, start_position, output_file):
    # Load fragment_data sheet from the xlsx file
    fragment_df = pd.read_excel(xlsx_file, sheet_name='fragment_data')
    
    # Load primer indexes from the csv file
    primer_df = pd.read_csv(csv_file)
    
    # Trim the CSV file to start at the specific position
    primer_df = primer_df.iloc[start_position:].reset_index(drop=True)
    
    # Create a dictionary for fast lookup of F_seq, R_seq_rc, and new_name based on name
    primer_map = {
        row['name']: (row['F_seq'], row['R_seq_rc'], row['new_name']) 
        for _, row in primer_df.iterrows()
    }
    
    def replace_oligo(row):
        oligo = row['oligo']
        # Extract the first 18 bases and last 18 bases from the oligo
        first_18 = oligo[:18]
        last_18 = oligo[-18:]
        
        replacement_done = False
        
        # Find matching primer set
        for name, (F_seq, R_seq_rc, new_name) in primer_map.items():
            if first_18 == F_seq and last_18 == R_seq_rc:
                # Check if it's already correctly formatted
                if oligo == (F_seq + oligo[18:-18] + R_seq_rc):
                    print(f"Oligo already correctly formatted for primer {name}. No change needed.")
                else:
                    # Replace the oligo's start and end with new sequences
                    new_oligo = F_seq + oligo[18:-18] + R_seq_rc
                    row['oligo'] = new_oligo
                    row['primer_index'] = new_name
                    replacement_done = True
                    print(f"Matched primer: {name}")
                    print(f"Original oligo: {oligo}")
                    print(f"New oligo: {new_oligo}")
                    print(f"Original primer_index: {row['primer_index']}")
                    print(f"New primer_index: {new_name}")
                break

        if not replacement_done:
            # Debugging output when no match is found
            print(f"No match found for oligo: {oligo}")
        
        return row
    
    # Apply the replacement across all rows
    fragment_df = fragment_df.apply(replace_oligo, axis=1)
    
    # Output to a new xlsx file
    fragment_df.to_excel(output_file, sheet_name='fragment_data', index=False)

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Reformat fragment data with updated primer indexes")
    parser.add_argument("xlsx_file", help="Input XLSX file with fragment data")
    parser.add_argument("csv_file", help="CSV file containing primer indexes")
    parser.add_argument("start_position", type=int, help="Starting position in the CSV primer index file")
    parser.add_argument("output_file", help="Output file name (e.g., input_reformatted.xlsx)")
    
    args = parser.parse_args()
    
    # Call the reformat function with command-line arguments
    reformat_fragments(args.xlsx_file, args.csv_file, args.start_position, args.output_file)
