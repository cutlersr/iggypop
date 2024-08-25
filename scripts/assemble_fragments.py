import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
from collections import defaultdict

"""
Iggypop Fragment Assembly Script

This script processes DNA oligo fragments from an iggypop FASTA file,
trims the ends to remove primers and cut site sequences, checks for matching 
overhangs between fragments, and assembles them into full sequences.

The main goals are to (i) confirm that a FASTA file containing gene fragment 
sequences can be properly assembled and (ii) ensure that the same indexing 
(barcode) primer set is not used more than once. This latter point is valuable 
to prevent errors when combining genes from different runs and users.

The script requires a standard nomenclature for the fragments:

>{accession}.{n}_FRAG_{m}_PRI_{set}

(n: Chisel number; m: Fragment number; set: primer set)

If you are analyzing iggypop output files, you don't need to change anything. 
It will ignore sequences that do not conform to that format.

Arguments:
--i            Input FASTA file
--o            Output FASTA file name
--n            Number of base pairs to trim from each end; for standard iggypop 
               is 25 (18 bp for primer + 7 bp for the BsmBI site). Use 26 for 
               two-step assemblies.
--cloning_ohs  External overhangs used for cloning (default: ['AATG', 'GCTT'])
--skip_library Skip assembling sequences with 'library' in the name. This is for 
               specialized projects that make libraries using iggypop.
"""


def read_fasta(input_file):
    """Reads sequences from a FASTA file and stores them in a dictionary."""
    sequences = {}
    for record in SeqIO.parse(input_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


def extract_base_gene_name(fragment_id):
    """
    Extracts the base gene name and primer information from the fragment ID 
    based on a predefined naming convention.
    """
    match = re.match(r"^(.*)_FRAG_.*_PRI_(.*\d+)", fragment_id)
    if match:
        base_gene_name = match.group(1)
        return base_gene_name, match.group(1), match.group(2), True
    return None, None, None, False


def process_sequences(sequences, n=25, skip_library=True):
    """
    Processes sequences by trimming bases, organizing them by gene name and 
    primer information, and optionally skipping "library" sequences.
    """
    processed_sequences = defaultdict(list)
    primer_pairs = defaultdict(set)
    all_primers = defaultdict(set)
    non_matching_count = 0
    skipped_library_count = 0
    library_primer_associations = defaultdict(set)

    for seq_id, seq in sequences.items():
        base_gene_name, gene_name, primer_info, is_match = extract_base_gene_name(seq_id)
        if is_match:
            if skip_library and "library" in seq_id.lower():
                skipped_library_count += 1
                library_primer_associations[primer_info].add(base_gene_name)
                continue

            processed_sequences[base_gene_name].append((gene_name, seq[n:-n]))
            primer_pairs[base_gene_name].add(primer_info)
            all_primers[primer_info].add(base_gene_name)
        else:
            non_matching_count += 1

    return (
        processed_sequences, primer_pairs, all_primers, non_matching_count, 
        skipped_library_count, library_primer_associations
    )


def assemble_sequences(processed_sequences):
    """
    Assembles sequences by concatenating trimmed fragments and checks for 
    mismatches between fragment ends.
    """
    assembled_sequences = {}
    overhangs = {}
    mismatches = defaultdict(list)

    for base_gene, fragments in processed_sequences.items():
        merged_fragments = defaultdict(list)
        for gene_name, fragment in fragments:
            merged_fragments[gene_name].append(fragment)

        for gene_name, gene_fragments in merged_fragments.items():
            assembled_seq = gene_fragments[0]
            gene_overhangs = [assembled_seq[:4]]  # First 4 bases
            for i in range(1, len(gene_fragments)):
                previous_fragment = gene_fragments[i - 1]
                current_fragment = gene_fragments[i]
                gene_overhangs.append(current_fragment[:4])  # Inner 4 bases
                if previous_fragment[-4:] != current_fragment[:4]:
                    mismatch_info = (
                        f"Warning: Mismatch between fragments {i-1} and {i} "
                        f"in {gene_name}: {previous_fragment[-4:]} != "
                        f"{current_fragment[:4]}"
                    )
                    mismatches[base_gene].append(mismatch_info)
                    print(mismatch_info)

                assembled_seq = assembled_seq[:-4] + current_fragment  # Assemble fragments
            gene_overhangs.append(assembled_seq[-4:])  # Last 4 bases of the assembled sequence
            if base_gene in assembled_sequences:
                assembled_sequences[base_gene] += assembled_seq
                overhangs[base_gene] += gene_overhangs[1:]
            else:
                assembled_sequences[base_gene] = assembled_seq
                overhangs[base_gene] = gene_overhangs

    return assembled_sequences, overhangs, mismatches


def write_fasta(
    assembled_sequences, primer_pairs, overhangs, all_primers, 
    library_primer_associations, cloning_ohs, output_file="assembled_seq.fasta"
):
    """
    Writes the assembled sequences to a FASTA file, including primer and overhang 
    information in the headers, and checks for primer reuse across multiple genes.
    """
    records = []
    for base_gene, seq in assembled_sequences.items():
        primers = " ".join(sorted(primer_pairs[base_gene]))  # Sort primers for consistency
        overhang_set = ", ".join(overhangs[base_gene])
        record_id = f"{base_gene} {primers} ohs:[{overhang_set}]"
        if not (seq.startswith(cloning_ohs[0]) and seq.endswith(cloning_ohs[1])):
            print(f"Warning: Cloning overhangs do not match for {record_id}")
        record = SeqRecord(Seq(seq), id=record_id, description="")
        records.append(record)

    # Check for primers used across multiple genes and print warnings
    problems_found = False
    for primer, genes in all_primers.items():
        if len(genes) > 1:
            print(f"Warning: Primer {primer} is used for multiple genes: {', '.join(genes)}")
            problems_found = True
        else:
            print(f"Primer {primer} is used for gene {', '.join(genes)}")

    # Check library primer associations
    for primer, genes in library_primer_associations.items():
        if len(genes) > 1:
            print(f"Warning: Library Primer {primer} is used for multiple genes: {', '.join(genes)}")
            problems_found = True
        else:
            print(f"Library Primer {primer} is used for gene {', '.join(genes)}")

    if not problems_found:
        print("All primers are used with only one gene/fragment_id.")

    SeqIO.write(records, output_file, "fasta")


def main():
    """Main function to parse command-line arguments and execute the script."""
    parser = argparse.ArgumentParser(description="Iggypop fragment assembler & checker")

    parser.add_argument("--i", help="Input FASTA file")

    parser.add_argument(
        "--n", type=int, default=25, 
        help=(
            "Number of base pairs to trim from each end; for standard iggypop this "
            "is 25 (18 bp for primer + 7 bp for the BsmBI site). Use 26 for two-step "
            "assemblies. Default: 25"
        )
    )
    parser.add_argument(
        "--o", default="assembled_seq.fasta", 
        help="Output FASTA file"
    )
    parser.add_argument(
        "--cloning_ohs", nargs=2, default=["AATG", "GCTT"], 
        help="Cloning overhangs (default: ['AATG', 'GCTT'])"
    )
    parser.add_argument(
        "--skip_library", action='store_true', default=True, 
        help=(
            "Skip assembling sequences with 'library' in the name. This is for "
            "specialized projects that make libraries using iggypop."
        )
    )

    args = parser.parse_args()

    sequences = read_fasta(args.i)
    (
        processed_sequences, primer_pairs, all_primers, non_matching_count, 
        skipped_library_count, library_primer_associations
    ) = process_sequences(sequences, args.n, args.skip_library)

    if not args.skip_library:
        print(f"Skipped {skipped_library_count} sequences containing 'library' in their name.")

    assembled_sequences, overhangs, mismatches = assemble_sequences(processed_sequences)
    write_fasta(
        assembled_sequences, primer_pairs, overhangs, all_primers, 
        library_primer_associations, args.cloning_ohs, args.o
    )

    if non_matching_count > 0:
        print(
            f"Warning: {non_matching_count} sequences did not match the expected format "
            "and were ignored."
        )


if __name__ == "__main__":
    main()
