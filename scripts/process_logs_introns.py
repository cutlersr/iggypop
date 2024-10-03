import re
import argparse
import pandas as pd
import os

def count_splices(donor_dict, acceptor_dict, min_distance=20):
    """
    Counts the number of splices, where a splice is defined as a donor
    occurring before an acceptor based on their integer positions and the distance
    between them is greater than the specified min_distance.
    """
    if not donor_dict or not acceptor_dict:
        return 0
    
    # Extract donor and acceptor positions
    donor_positions = sorted(map(int, donor_dict.keys()))
    acceptor_positions = sorted(map(int, acceptor_dict.keys()))

    # Count splices where a donor occurs before an acceptor with a distance greater than min_distance
    splices = 0
    for donor in donor_positions:
        for acceptor in acceptor_positions:
            if donor < acceptor and (acceptor - donor) > min_distance:
                splices += 1

    return splices

def parse_sequence_name(sequence_name):
    """
    Parses the sequence name into three fields: accession, repeat, and chisel.
    If the chisel component is missing, it defaults to 1.
    """
    sequence_name = sequence_name.lstrip('>')  # Remove the leading ">" character
    match = re.match(r'(?P<accession>.+)\.(?P<repeat>\d+)_chisel(?:_(?P<chisel>\d+))?', sequence_name)
    
    if match:
        accession = match.group('accession')
        repeat = match.group('repeat')
        chisel = match.group('chisel') if match.group('chisel') else '1'
        return accession, repeat, chisel
    else:
        return sequence_name, None, '1'  # In case the pattern doesn't match, default chisel to 1

def extract_sequence(content, sequence_name):
    """
    Extracts the multiline sequence that follows the sequence name and stops before the keyword "Chiseled".
    """
    # Regular expression to capture the sequence block
    sequence_pattern = re.compile(
        rf'>{re.escape(sequence_name)}\n(?P<sequence>[ATGC\n]+)\nChiseled', 
        re.DOTALL
    )
    
    # Search for the sequence block
    match = sequence_pattern.search(content)
    if match:
        sequence = match.group('sequence').replace('\n', '')  # Remove newlines
        return sequence
    return None

def extract_blocks(log_file, min_distance=20):
    with open(log_file, 'r') as file:
        content = file.read()

    # Regular expression for sequence names with chisel
    sequence_pattern = re.compile(r'>.*_chisel(?:_\d)?')

    # Regular expression for the blocks of interest
    block_pattern = re.compile(
        r'(scanning chiseled sequence for cryptic introns.*?)(?:possible cryptic intron found|no intron found)', 
        re.DOTALL
    )

    # Regular expression for intron status
    intron_pattern = re.compile(r'intron:\s(True|False)')

    # Regular expression for donor block
    donor_pattern = re.compile(r'dico_donor:\s({.*?})', re.DOTALL)

    # Regular expression for acceptor block
    acceptor_pattern = re.compile(r'dico_acceptor:\s({.*?})', re.DOTALL)

    # Find all sequence names
    sequence_names = sequence_pattern.findall(content)

    # Find all blocks
    blocks = block_pattern.findall(content)

    # Initialize a list to store the extracted data for the DataFrame
    data = []

    # Iterate through each block and extract relevant info
    for i, seq in enumerate(sequence_names):
        accession, repeat, chisel = parse_sequence_name(seq)
        sequence_content = extract_sequence(content, seq.lstrip('>'))  # Extract the sequence

        if i < len(blocks):
            block = blocks[i]
            
            # Check intron status
            intron_match = intron_pattern.search(block)
            intron_status = intron_match.group(1) if intron_match else "False"
            
            # Extract donor block if present
            donor_match = donor_pattern.search(block)
            donor_block = donor_match.group(1) if donor_match else "n/a"
            
            # Extract acceptor block if present
            acceptor_match = acceptor_pattern.search(block)
            acceptor_block = acceptor_match.group(1) if acceptor_match else "n/a"

            # Convert donor and acceptor blocks to dicts if they're found
            if donor_block != "n/a":
                donor_dict = eval(donor_block)
            else:
                donor_dict = None
            
            if acceptor_block != "n/a":
                acceptor_dict = eval(acceptor_block)
            else:
                acceptor_dict = None

            # Count the number of entries in the donor and acceptor dictionaries
            donor_count = len(donor_dict) if donor_dict else 0
            acceptor_count = len(acceptor_dict) if acceptor_dict else 0

            # Count valid splices with the condition of donor < acceptor and distance > min_distance
            splices = count_splices(donor_dict, acceptor_dict, min_distance)
            
            # Append the extracted data as a row to the list
            data.append({
                'accession': accession,
                'repeat': int(repeat),
                'chisel': int(chisel),
                'intron': intron_status == "True",  # Store as boolean
                'donor_block': donor_dict,
                'acceptor_block': acceptor_dict,
                'donors': donor_count,
                'acceptors': acceptor_count,
                'splices': splices,
                'sequence': sequence_content  # Add the extracted sequence
            })
        else:
            # Append default values for sequences without a corresponding block
            data.append({
                'accession': accession,
                'repeat': int(repeat),
                'chisel': int(chisel),
                'intron': False,  # No intron by default
                'donor_block': None,
                'acceptor_block': None,
                'donors': 0,
                'acceptors': 0,
                'splices': 0,
                'sequence': sequence_content  # Add the extracted sequence
            })

    # Create a DataFrame from the collected data
    df = pd.DataFrame(data)

    # Write the DataFrame to a CSV file
    output_csv = os.path.join(os.path.dirname(log_file), 'intron_data.csv')
    df.to_csv(output_csv, index=False)

    print(f"Data has been written to {output_csv}")
    print(f"Repeats processed: {df['repeat'].nunique()}")
    print(f"Total sequences with splices: {(df['splices'] > 0).sum()}")
    print(f"Total sequences without splices: {(df['splices'] == 0).sum()}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract chiseled sequence blocks from log file.")
    parser.add_argument('log_file', type=str, help="Path to the log file")
    parser.add_argument('--min_distance', type=int, default=20, help="Minimum distance between donor and acceptor")
    args = parser.parse_args()

    extract_blocks(args.log_file, min_distance=args.min_distance)
