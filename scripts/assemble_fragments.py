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

For one-step assemblies:
  >{accession}.{n}_FRAG_{m}_PRI_{set}

For two-step assemblies (which require a trim of 26 in step 1), the sequence 
name should include a pattern like "L0_<digit>" (e.g.:
  >STARBURST_hispS_npgA_H3H_Luz_CPH_v1_npgA_removed.1_L0_3_1_PRI_set234)

In two-step assemblies the base gene is recognized as everything up to (and including) 
the "L0_<digit>" (e.g. STARBURST_hispS_npgA_H3H_Luz_CPH_v1_npgA_removed.1_L0_3) in step 1.
Then, in a second assembly step, all first-step assemblies from the same two-step group 
are reassembled (with a trim of 25 bp) to yield the final gene, e.g.:
  STARBURST_hispS_npgA_H3H_Luz_CPH_v1_npgA_removed.1

The final output FASTA file will include:
  - One-step assemblies (for fragments not in two-step groups)
  - For two-step groups, both the first-step assembly (prefixed with "STEP1:")
    and the final second-step assembly (prefixed with "STEP2:").

Arguments:
--i            Input FASTA file
--o            Output FASTA file name (final assemblies)
--n            Number of bp to trim from each end for 1-step assemblies; 
               default is 25; sequences with L0_<digit> are trimmed with 26 bp in step 1.
--cloning_ohs  External overhangs used for cloning (default: ['AATG', 'GCTT'])
--skip_library Skip assembling sequences with 'library' in the name.
"""

def read_fasta(input_file):
    """Reads sequences from a FASTA file and returns a dict {record.id: record.seq}."""
    sequences = {}
    for record in SeqIO.parse(input_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def extract_base_gene_name(fragment_id):
    """
    Extracts the base gene name and primer info from the fragment ID.
    
    For two-step assemblies (IDs containing "L0_<digit>"), captures everything 
    up to the "_L0_<digit>" portion. For one-step sequences, uses the standard pattern.
    """
    if re.search(r"L0_\d", fragment_id):
        # Example: STARBURST..._L0_3_1_PRI_set234 -> returns key "STARBURST..._L0_3"
        match = re.match(r"^(.*_L0_\d+)_\d+_PRI_(.*\d+)", fragment_id)
    else:
        match = re.match(r"^(.*)_FRAG_.*_PRI_(.*\d+)", fragment_id)
    if match:
        base_gene_name = match.group(1)
        primer_info = match.group(2)
        return base_gene_name, base_gene_name, primer_info, True
    return None, None, None, False

def process_sequences(sequences, n=25, skip_library=True):
    """
    Processes sequences by trimming and grouping by gene name.
    
    For two-step fragments (ID contains L0_<digit>), uses a trim of 26 bp (step 1).
    Also collects oligo statistics.
    """
    processed_sequences = defaultdict(list)
    primer_pairs = defaultdict(set)
    all_primers = defaultdict(set)
    library_primer_associations = defaultdict(set)
    primer_types = {}  # Map base_gene to type ("FRAG" or "LIBRARY")
    
    total_oligo_count = 0
    gene_oligo_count = 0
    library_oligo_count = 0
    gene_oligo_length_sum = 0
    library_oligo_length_sum = 0
    non_matching_count = 0
    skipped_library_count = 0
    non_matching_ids = []
    
    for seq_id, seq in sequences.items():
        base_gene_name, gene_name, primer_info, is_match = extract_base_gene_name(seq_id)
        if is_match:
            total_oligo_count += 1
            current_type = "LIBRARY" if "library" in seq_id.lower() else "FRAG"
            primer_types[base_gene_name] = current_type
            trim_val = 26 if re.search(r"L0_\d", seq_id) else n
            trimmed_seq = seq[trim_val:-trim_val]
            if current_type == "LIBRARY":
                library_oligo_count += 1
                library_oligo_length_sum += len(trimmed_seq)
            else:
                gene_oligo_count += 1
                gene_oligo_length_sum += len(trimmed_seq)
            if skip_library and current_type == "LIBRARY":
                skipped_library_count += 1
                library_primer_associations[primer_info].add(base_gene_name)
                continue
            processed_sequences[base_gene_name].append((gene_name, trimmed_seq))
            primer_pairs[base_gene_name].add(primer_info)
            all_primers[primer_info].add(base_gene_name)
        else:
            non_matching_count += 1
            non_matching_ids.append(seq_id)
    
    return (processed_sequences, primer_pairs, all_primers, non_matching_count,
            skipped_library_count, library_primer_associations, non_matching_ids,
            primer_types, total_oligo_count, gene_oligo_count, library_oligo_count,
            gene_oligo_length_sum, library_oligo_length_sum)

def assemble_sequences(processed_sequences):
    """
    Assembles sequences (first-step) by concatenating trimmed fragments.
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
            gene_overhangs = [assembled_seq[:4]]
            for i in range(1, len(gene_fragments)):
                prev_frag = gene_fragments[i - 1]
                curr_frag = gene_fragments[i]
                gene_overhangs.append(curr_frag[:4])
                if prev_frag[-4:] != curr_frag[:4]:
                    msg = (f"Warning: Mismatch between fragments {i-1} and {i} "
                           f"in {gene_name}: {prev_frag[-4:]} != {curr_frag[:4]}")
                    mismatches[base_gene].append(msg)
                    print(msg)
                assembled_seq = assembled_seq[:-4] + curr_frag
            gene_overhangs.append(assembled_seq[-4:])
            assembled_sequences[base_gene] = assembled_seq
            overhangs[base_gene] = gene_overhangs
    return assembled_sequences, overhangs, mismatches

def simulate_second_step(first_step, n_second=25, cloning_ohs=["AATG", "GCTT"]):
    """
    Simulates the second assembly step for two-step fragments.
    
    For each key in first_step that contains "_L0_", group all such keys by the final gene name
    (obtained by stripping the "_L0_*" suffix). Then, sort the fragments (using the numeric value
    following "L0_" if available), concatenate them using the overlap logic, trim n_second bp from
    each end, and enforce that the final sequence starts with cloning_ohs[0] and ends with cloning_ohs[1].
    """
    groups = defaultdict(list)
    final_assemblies = {}
    
    # Partition keys: if key contains "_L0_", group by final gene name.
    for key, seq in first_step.items():
        if "_L0_" in key:
            group_key = re.sub(r"_L0_.*", "", key)
            m = re.search(r"_L0_(\d+)", key)
            num = int(m.group(1)) if m else 0
            groups[group_key].append((num, seq))
    # Assemble each group.
    for group_key, frag_list in groups.items():
        frag_list.sort(key=lambda x: x[0])
        assembled = frag_list[0][1]
        for _, frag in frag_list[1:]:
            if assembled[-4:] != frag[:4]:
                print(f"Warning (step2): Overlap mismatch in group {group_key}: {assembled[-4:]} != {frag[:4]}")
            assembled = assembled[:-4] + frag
        if len(assembled) > 2 * n_second:
            final_seq = assembled[n_second:-n_second]
        else:
            final_seq = assembled
        if not final_seq.startswith(cloning_ohs[0]):
            final_seq = cloning_ohs[0] + final_seq
        if not final_seq.endswith(cloning_ohs[1]):
            final_seq = final_seq + cloning_ohs[1]
        final_assemblies[group_key] = final_seq
    return final_assemblies

def write_fasta(assemblies, primer_pairs, overhangs, all_primers,
                library_primer_associations, cloning_ohs, output_file):
    """
    Writes the assembled sequences to a FASTA file.
    """
    records = []
    for gene, seq in assemblies.items():
        primers = " ".join(sorted(primer_pairs.get(gene, [])))
        ohs = ", ".join(overhangs.get(gene, []))
        record_id = f"{gene} {primers} ohs:[{ohs}]"
        if not (seq.startswith(cloning_ohs[0]) and seq.endswith(cloning_ohs[1])):
            print(f"Warning: Cloning overhangs do not match for {record_id}")
        record = SeqRecord(Seq(seq), id=record_id, description="")
        records.append(record)
    for primer, genes in all_primers.items():
        if len(genes) > 1:
            print(f"Warning: Primer {primer} is used for multiple genes: {', '.join(genes)}")
    for primer, genes in library_primer_associations.items():
        if len(genes) > 1:
            print(f"Warning: Library Primer {primer} is used for multiple genes: {', '.join(genes)}")
    SeqIO.write(records, output_file, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Iggypop fragment assembler & checker")
    parser.add_argument("--i", help="Input FASTA file", required=True)
    parser.add_argument("--n", type=int, default=25,
                        help=("Number of bp to trim from each end for 1-step assemblies; "
                              "default is 25; sequences with L0_<digit> are trimmed with 26 bp in step 1."))
    parser.add_argument("--o", default="assembled_seq.fasta", help="Output FASTA file (final assemblies)")
    parser.add_argument("--cloning_ohs", nargs=2, default=["AATG", "GCTT"],
                        help="Cloning overhangs (default: ['AATG', 'GCTT'])")
    parser.add_argument("--skip_library", action='store_true', default=True,
                        help="Skip assembling sequences with 'library' in their name.")
    args = parser.parse_args()

    sequences = read_fasta(args.i)
    (processed_sequences, primer_pairs, all_primers, non_matching_count,
     skipped_library_count, library_primer_associations, non_matching_ids,
     primer_types, total_oligo_count, gene_oligo_count, library_oligo_count,
     gene_oligo_length_sum, library_oligo_length_sum) = process_sequences(sequences, args.n, args.skip_library)

    if non_matching_count > 0:
        print(f"Warning: {non_matching_count} sequences did not match the expected format and were ignored.")
        for seq_id in non_matching_ids:
            print(f"  - {seq_id}")
    if not args.skip_library and skipped_library_count:
        print(f"Skipped {skipped_library_count} sequences containing 'library' in their name.")

    # First-step assembly.
    first_step_assemblies, first_overhangs, mismatches = assemble_sequences(processed_sequences)
    
    # Separate one-step and two-step first-step assemblies.
    one_step_assemblies = {k: v for k, v in first_step_assemblies.items() if "_L0_" not in k}
    twostep_first = {k: v for k, v in first_step_assemblies.items() if "_L0_" in k}
    
    # Simulate second-step assembly for two-step groups.
    twostep_final = simulate_second_step(first_step_assemblies, n_second=25, cloning_ohs=args.cloning_ohs)
    
    # Create combined dictionary with labels.
    final_dict = {}
    final_primer_pairs = {}
    final_overhangs = {}
    
    # One-step assemblies: keys unchanged.
    for k, seq in one_step_assemblies.items():
        final_dict[k] = seq
        final_primer_pairs[k] = primer_pairs[k]
        final_overhangs[k] = first_overhangs[k]
    
    # Two-step first-step assemblies: add prefix "STEP1:".
    for k, seq in twostep_first.items():
        new_key = "STEP1:" + k
        final_dict[new_key] = seq
        final_primer_pairs[new_key] = primer_pairs[k]
        final_overhangs[new_key] = first_overhangs[k]
    
    # Two-step final assemblies: add prefix "STEP2:".
    for k, seq in twostep_final.items():
        new_key = "STEP2:" + k
        final_dict[new_key] = seq
        final_primer_pairs[new_key] = set()  # Not available for final step.
        final_overhangs[new_key] = []       # Not available for final step.
    
    # Write all combined assemblies to one output FASTA file.
    write_fasta(final_dict, final_primer_pairs, final_overhangs,
                all_primers, library_primer_associations, args.cloning_ohs, args.o)
    
    # Statistics.
    total_final = len(final_dict)
    total_bases = sum(len(seq) for seq in final_dict.values())
    avg_final_length = total_bases / total_final if total_final > 0 else 0
    overall_oligo_count = total_oligo_count
    overall_length_sum = gene_oligo_length_sum + library_oligo_length_sum
    avg_oligo_length = overall_length_sum / overall_oligo_count if overall_oligo_count > 0 else 0
    
    # Two-step stats.
    two_step_keys = [k for k in primer_types if "_L0_" in k]
    step1_groups = len(two_step_keys)
    final_twostep = sum(1 for k in final_dict if k.startswith("STEP2:"))
    
    print("\n--- Assembly Statistics ---")
    print(f"Total assembled genes (final): {total_final}")
    print(f"Average assembled sequence length: {avg_final_length:.2f} bp")
    print(f"Total number of bases synthesized: {total_bases}")
    print(f"Total number of oligonucleotides processed: {overall_oligo_count}")
    print(f"Total oligos for genes: {gene_oligo_count}")
    print(f"Total oligos for libraries: {library_oligo_count}")
    print(f"Average oligo length: {avg_oligo_length:.2f} bp")
    print("\n--- Two-step Assembly Stats ---")
    print(f"Step 1 two-step groups: {step1_groups}")
    print(f"Step 2 assembled two-step genes: {final_twostep}")

if __name__ == "__main__":
    main()
