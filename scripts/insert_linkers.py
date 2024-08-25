import argparse
import csv

"""
 This script generates new iggypop fragments by taking an existing fragment and inserting 
 sequences from a fasta file into regions marked with []. It additionally picke new barcodes
 The new sequences and used barcodes are saved to output files. We have used this to design
 fragment libraries with differing linker variants at specified points. Each linker library
 designed can be amplified with a unique barcode and then used together with other iggypop
 fragments to synthesize a library of variants.

"""

def read_fasta(filename):
    """Reads a FASTA file and returns a dictionary with sequence names as keys and sequences as values."""
    sequences = {}
    with open(filename, 'r') as file:
        name = None
        seq = []
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if name:
                    sequences[name] = ''.join(seq)
                name = line[1:]
                seq = []
            else:
                seq.append(line)
        if name:
            sequences[name] = ''.join(seq)
    return sequences

def read_barcodes(filename):
    """Reads a CSV file with barcodes and returns a list of tuples (name, F_seq, R_seq_rc)."""
    barcodes = []
    with open(filename, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            barcodes.append((row['name'], row['F_seq'], row['R_seq_rc']))
    return barcodes

def generate_new_sequences(input_sequences, linkers, barcodes, start_index):
    """Generates new sequences by replacing the region in [] with each linker and adding barcodes."""
    new_sequences = {}
    used_barcodes = []
    barcode_index = start_index
    for input_name, input_seq in input_sequences.items():
        if '[' in input_seq and ']' in input_seq:
            pre_region, post_region = input_seq.split('[', 1)
            region_to_replace, post_region = post_region.split(']', 1)
            barcode = barcodes[barcode_index % len(barcodes)]
            barcode_name, F_seq, R_seq_rc = barcode
            for linker_name, linker_seq in linkers.items():
                new_name = f"{input_name}_{linker_name}_{barcode_name}"
                new_seq = f"{F_seq}{pre_region}{linker_seq}{post_region}{R_seq_rc}"
                new_sequences[new_name] = new_seq
            used_barcodes.append((barcode_index, barcode_name))
            barcode_index += 1
        else:
            # If no region to replace, just add the original sequence with barcode
            barcode = barcodes[barcode_index % len(barcodes)]
            barcode_name, F_seq, R_seq_rc = barcode
            new_name = f"{input_name}_{barcode_name}"
            new_seq = f"{F_seq}{input_seq}{R_seq_rc}"
            new_sequences[new_name] = new_seq
            used_barcodes.append((barcode_index, barcode_name))
            barcode_index += 1
    return new_sequences, used_barcodes

def write_fasta(sequences, filename):
    """Writes sequences to a FASTA file."""
    with open(filename, 'w') as file:
        for name, seq in sequences.items():
            file.write(f">{name}\n")
            for i in range(0, len(seq), 80):  # Wrap lines at 80 characters
                file.write(f"{seq[i:i+80]}\n")

def write_used_barcodes(used_barcodes, filename):
    """Writes used barcodes and their indices to a TXT file."""
    with open(filename, 'w') as file:
        for index, name in used_barcodes:
            file.write(f"Index: {index}, Barcode: {name}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Generate sequences by inserting linkers into regions marked with [].")
    parser.add_argument(
        '--i', required=True, 
        help="Input FASTA file with sequences containing regions in []."
    )
    parser.add_argument(
        '--linkers', required=True, 
        help="FASTA file with linkers."
    )
    parser.add_argument(
        '--barcode', required=True, 
        help="CSV file with barcodes."
    )
    parser.add_argument(
        '--index', type=int, default=0, 
        help="Row number of the barcode CSV to start adding the barcodes from "
             "(0-based index)."
    )
    args = parser.parse_args()

    input_sequences = read_fasta(args.i)
    linkers = read_fasta(args.linkers)
    barcodes = read_barcodes(args.barcode)

    new_sequences, used_barcodes = generate_new_sequences(
        input_sequences, linkers, barcodes, args.index)

    output_file = 'output.fasta'
    barcodes_file = 'used_barcodes.txt'
    write_fasta(new_sequences, output_file)
    write_used_barcodes(used_barcodes, barcodes_file)
    print(f"Generated sequences written to {output_file}")
    print(f"Used barcodes written to {barcodes_file}")
    for index, name in used_barcodes:
        print(f"Index: {index}, Barcode: {name}")

if __name__ == "__main__":
    main()
