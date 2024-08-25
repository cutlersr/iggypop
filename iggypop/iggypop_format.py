import argparse
import yaml
from pop_helpers import annotate_genbank, print_organisms_from_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="pop_format.py adds pop constraint and optimization instructions to a GenBank record using a yaml configuration file.")
    parser.add_argument("--i", default="test.gb", help="Input GenBank file/path")
    parser.add_argument("--o", default="modified_genbank.gb", help="output file name/path")
    parser.add_argument("--pct", default=20, type=int, help="setting for EnforceChanges to use with hybrid codon opt; default=20")
    parser.add_argument("--id", default=False, help="n of copies to make; useful when you are tuning by tweaking")
    parser.add_argument("--yml", default="yaml/gb_mcu.yml", help="Path to the input YAML file")
    parser.add_argument("--taxIDs", action="store_true", help="Print the list of common organisms and their NCBI taxIDs")
    
    args = parser.parse_args()
   
    if args.taxIDs:
        print_organisms_from_file('data/organisms.txt')
        exit(0)  # Exit after printing the organisms

    annotate_genbank(args.i, args.o, args.yml, args.id, args.pct)
 