"""
This program is designed to create barcoded fragments of target coding
sequences that can be amplified in a single PCR reaction and then assembled
using Golden Gate cloning. Input sequences are optionally optimized using
DNAchisel. Depending on the parameters set, you can domesticate sequences
to remove restriction enzymes, perform codon optimization, prevent hairpins,
and apply various other optimization tools offered by DNAchisel prior to
generating the fragmeted/barcoded sequences.

DNAchisel's three codon optimization methods can be used (`match_codon_usage`,
`use_best_codon`, `harmonize_rca`). Additionally, a "hybrid" mode is available
that uses `use_best_codon` with a constraint that the sequence be `--pct`
different from the input; this leads to higher CAI values than
`match_codon_usage` and sequences with greater diversity than `use_best_codon`.

An experimental feature (`--deintronize`) uses Spliceator in combination with
DNAchisel to remove potential cryptic introns and splice donor/acceptor
sites. If you want to skip chiseling or sequence optimization, use
`--mode no_mods`.

Once optimized, sequences are appended with defined 5' and 3' ends to create_
overhangs and other elements needed for Golden Gate cloning. The default
settings add sequences generate a MoClo-compatible ORF, resulting in
final clones flanked by BsaI sites and AATG/GCTT overhangs (after assembly), 
thus suitable to build Level 1 constructs in MoClo-compatible systems.

The final target sequences are fragmented into segments with high-fidelity
overhangs for efficient reassembly via Golden Gate cloning. The maximum
fragment size is set by the `--segment_length` flag, with a default value
of 200, which is appropriate for synthesizing 250 bp oligos (reserving 50 bp
for enzyme sites and barcode primers). Fragmentation is done with
GoldenHinges, which limits overhang selection to predefined sets. To
ensure high-fidelity assemblies, sets of precomputed high-fidelity overhangs
are used as constraints. Multiple overhang sets are considered (`--n_tries`),
with the highest fidelity set reported.

Assemblies for sequences larger than a few kilobases can be assembled in
a two-step process (`--two_step on` flag), aiming for ~1 kb first-step
assemblies. This approach increases the frequency of error-free clones
required in the second step.

By default, a single barcoding primer set is appended to each fragment,
enabling a single-tube reaction for reassembly. For more complex assemblies,
multiple barcodes can be used per gene, with the `--max_fragments` flag
determining the number of barcodes. Barcodes are added from a file in the
`data/` folder; you can use a custom set of barcodes if required with the
`--index_primers` flag.

Additional functionality includes the `--repeats` flag, which allows
generating multiple versions of a single input sequence, useful for testing
multiple optimized sequences. When combined with `--tweak_n`, it generates
a large set of candidates, from which a maximally diverse subset can be
selected for further testing.

Parameters can be set via command line or more conveniently with a YAML
file. Command-line arguments override YAML settings.

**Command line options**:
- `--i` (str): Input FASTA file path for sequences to optimize and fragment.
- `--o` (str): Output file stem, saved in the `out/` folder.
- `--yml` (str): YAML file for run parameters (default: `yaml/moclo_cds.yml`).
  Command-line arguments override YAML settings.
- `--mode` (str): Operation mode:
  - `chisel`: Modifies the input sequence (default).
  - `no_mods`: Only hinges and barcodes fragments.
  - `no_hinge`: Chisels sequences without hinging/barcoding.
- `--reports` (bool): Enables DNAchisel reports.

**Codon Optimization**:
- `--codon_opt` (str): Codon optimization method (`use_best_codon`,
  `match_codon_usage`, `harmonize_rca`, `hybrid`, `none`).
- `--pct` (float): Target sequence divergence for hybrid optimization
  (default: 20%).
- `--species` (str): Codon table for optimization (default: `arabidopsis`).
  Supports NCBI taxIDs.
- `--codon_tbl` (str): Codon table source (`cocoputs`, `kazusa`).
- `--original_species` (str): Source species for codon harmonization.
- `--deintronize` (str): Enables experimental deintronization mode.

**Repeated Design**:
- `--repeats` (int): Number of chiseled sequences to create per input.
- `--tweak_n` (int): Number of tweaked sequences for diversity.
- `--tweak_cai` (float): Minimum CAI value for tweaking (default: 0.8).

**Assembly**:
- `--two_step` (str): Enables two-step assemblies.
- `--max_fragments` (int): Maximum fragments per PCR (default: 6).
- `--segment_length` (int): Maximum segment length (default: 200 bp).
- `--ext_overhangs` (list): External overhangs for cloning, excluded from
  internal junctions.
- `--base_5p_end`, `--base_3p_end` (str): Sequences appended to the 5' and
  3' ends of the chiseled CDS.
- `--pcr_5p_cut`, `--pcr_3p_cut` (str): Sequences added to oligos for Golden
  Gate cloning.
- `--primer_index` (int): Starting point for adding primers from the index
  set.
- `--n_tries` (int): Number of overhang sets to consider (default: 50).
- `--radius` (int): Distance from ideal cut sites for selecting overhangs
  (default: 8).

**Miscellaneous**:
- `--seed` (int): Seed for random number generation.
- `--index_primers` (str): Path to index primers file (`data/10K_primers_renamed.csv`).
- `--fidelity_data` (str): Path to OH fidelity data (`data/FileS03_T4_18h_25C.xlsx`).
- `--ohsets` (str): High-fidelity overhang sets for assembly.
- `--taxIDs` (bool): Prints a list of common organisms and their NCBI taxIDs.
"""

import os
import sys
import datetime
import subprocess
import numpy
import random
import math
from headers import *
from initialization import *
from chisel_hinge import *
import pandas as pd
import yaml
import warnings
from pop_helpers import *
from python_codon_tables import get_codons_table
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Ignore warnings
warnings.simplefilter(action='ignore', category=Warning)

def hinge(chiseled_sequence, seq_id, segment_length):
    """
    Process a sequence to find the best hinge solution.

    Parameters:
    - chiseled_sequence: Modified sequence after codon optimization.
    - seq_id: Identifier for the sequence.
    - segment_length: Desired length of sequence segments.

    Returns:
    - A dictionary containing the best solutions for the sequence.
    """
    try:
        best_solution, seq_sets, fidelity_scores, \
            counter_value_for_best_solution = find_cut_solution(
                chiseled_sequence, overhang_sets, radius, ext_overhangs,
                segment_length, n_tries, potapov_data
            )

        best_overall_solution = max(fidelity_scores.values())
        log_and_print(
            f"\nOH set fidelity: {best_overall_solution:.3f}\n", log_file, quiet
        )

        max_fidelity_index = max(fidelity_scores, key=fidelity_scores.get)
        solution_data = seq_sets[max_fidelity_index]
        solutions = solution_data['solutions']
        counter_value = solution_data['counter']
        fidelity = fidelity_scores[max_fidelity_index]

        df = create_fragments_df(
            seq_id, solutions, original_sequence, chiseled_sequence,
            fidelity, base_5p_end, base_3p_end, original_cai, chiseled_cai
        )

        best_solutions_dict[seq_id] = {
            'dataframe': df,
            'fidelity': fidelity,
            'Tries_for_hinges': counter_value
        }

        return best_solutions_dict

    except Exception as e:
        log_and_print(
            f"{seq_id} FAILED hinge.\nException occurred: {e}\n", log_file
        )
        return None


if __name__ == "__main__":

    try:

########     
        tag, ofile, log_file_path, updated_defaults = initialize()
        globals().update(vars(updated_defaults))
        log_file = open(log_file_path, "a")

        log_file.write(f"Log file for: {tag}\n")
        log_file.write(
            f"Start time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        )
        log_file.write(f"Command line: {' '.join(sys.argv)}\n")

        print_header_cds(tag)

        if check_ext_overhangs(ext_overhangs):
            sys.exit()

        # Load Potapov data for calculating overhang fidelities
        try:
            potapov_data = pd.read_excel(fidelity_data)
        except Exception as e:
            log_and_exit(
                f"FAILED to load fidelity data.\nException occurred: {e}\n",
                log_file, exit_code=2
            )

        # Load high-fidelity overhang sets
        try:
            overhang_sets = get_overhang_sets(ohsets, ext_overhangs)
        except Exception as e:
            log_and_exit(
                f"FAILED to get overhang sets.\nException occurred: {e}\n",
                log_file, exit_code=3
            )

        # Read input sequences
        sequences = read_fasta(i)
        if not sequences:
            log_and_exit(
                "No sequences found in the input file.", log_file, exit_code=4
            )

        # Convert species to NCBI taxID for dnachisel
        if species == "arabidopsis":
            species = '3702'

        # Load codon table for specified species using NCBI taxID
        try:
            if (codon_tbl=="cocoputs"):

                if species == "arabidopsis":
                    species = 'a_thaliana'
                    print("using AT")

                codon_table, species_name = calculate_codon_frequencies(
                    file_path='data/cleaned_coco.tsv',
                    species_identifier=species
                )

            else:
                codon_table = get_codons_table(species)

            log_and_print(f"using {codon_tbl} for your codon table", log_file)
            log_and_print(f"using {codon_opt} for your codon optimization", log_file)


        except Exception as e:
            log_and_exit(
                "FAILED codon table lookup; this may be because you've used "
                "the wrong species name, taxID, do not have a network connection,"
                "or http://www.kazusa.or.jp/codon is down.\n"
                "Use the --codon_tbl cocoputs flag to load tables locally or, or use" 
                "`--codon_opt none' to skip codon lookup and optimization.",
                log_file,
                exit_code=5
            )

        if (deintronize=="on"):
            log_and_print("Running deintronize mode\n", log_file)
            try:
                from intron_stuff import *
            except Exception as e:
                log_and_print(
                    f"Error loading intron_stuff.\nException occurred: {e}\n"
                    "Iggypop will continue without deintronization enabled",
                    log_file
                )

        if taxIDs:
            print_taxIDs()
            sys.exit(0)

#####

        if tweak_n and mode != "no_hinge":
            print("Tweaking only works in no_hinge mode; run with '--mode no_hinge'.")            
            sys.exit()

        if tweak_n and repeats < tweak_n:
            print(f"Tweaking requires the use of repeats; run with \
                '--repeats {2 * tweak_n}' (or greater).")
            sys.exit()


        df_all = pd.DataFrame()
        two_step_df = pd.DataFrame()

#       Main sequence modification and fragmentation loop

        for i, seq_dict in enumerate(sequences):
            original_sequence = seq_dict['sequence'].upper()
            seq_id = seq_dict['id']

            seed = seed if seed else random.randint(0, 2**32 - 1)
            numpy.random.seed(seed)


            for repeat in range(repeats):
                accession = seq_id
                current_seq_id = f"{seq_id}.{repeat + 1}"
                best_solutions_dict = {}
                best_solution = None
                best_fidelity = 0
                intron_constraints = []
                avoid_changes_constraints = []
                loc_constraints_left = []
                loc_constraints_right = []
                original_cai = None
                chiseled_cai = None
                reports = reports

                log_and_print(f'seed: {seed}', log_file, quiet)

                log_and_print(f'\nProcessing {accession}', log_file)
                log_and_print(f"Seed value for this run: {seed}", log_file, quiet)
                log_and_print('*' * 80, log_file)

                if mode != 'no_mods':
                    if not check_orf(
                        original_sequence, seq_id, allowed_chars,
                        segment_length, log_file
                    ):
                        break

                    if codon_opt == 'harmonize_rca' and original_species == 'none':
                        log_and_exit(
                            "harmonize_rca requires that you specify the "
                            "CDS's origin species", log_file, exit_code=7
                        )

                    last_try = 0
                    max_attempts = 5
                    donor_seqs = []
                    acceptor_seqs = []
                    donor_locs_left = []
                    acceptor_locs_left = []
                    donor_locs_right = []
                    acceptor_locs_right = []
                    left_bounds = []
                    right_bounds = []

                    while last_try < max_attempts:
                        try:
                            if last_try == 0:
                                log_and_print(
                                    f">{current_seq_id}\n{original_sequence}",
                                    log_file, quiet
                                )

                            results_path = (
                                f'out/{tag}/{current_seq_id}_chisel_'
                                f'{last_try + 1}'
                            )

                            chiseled_sequence = chisel(
                                original_sequence, codon_opt, codon_table,
                                original_species, intron_constraints,
                                loc_constraints_left, loc_constraints_right,
                                left_bounds, right_bounds, log_file,
                                reports, pct, last_try, repeat, quiet, yml, results_path
                            )
                            time.sleep(0.2)

                            if last_try == 0:
                                log_and_print(
                                    f">{current_seq_id}_chisel\n"
                                    f"{chiseled_sequence}",
                                    log_file, quiet
                                )

                            if last_try > 0:
                                log_and_print(
                                    f">{current_seq_id}_chisel_{last_try + 1}\n"
                                    f"{chiseled_sequence}", log_file, quiet
                                )

                            if codon_opt != "none":
                                original_cai = calculate_cai(
                                    original_sequence, codon_table
                                )
                                chiseled_cai = calculate_cai(
                                    chiseled_sequence, codon_table
                                )

                                if codon_tbl == "kazusa":
                                    species_name = (
                                        "arabidopsis" if species == '3702'
                                        else str(species)
                                    )

                                log_and_print(
                                    f"\nChiseled with {codon_opt} optimization "
                                    f"and the {species_name} codon table.",
                                    log_file, quiet
                                )
                                log_and_print(
                                    f"\tOriginal CAI: {original_cai:.3f}", log_file, quiet
                                )
                                log_and_print(
                                    f"\tChiseled CAI: {chiseled_cai:.3f}\n",
                                    log_file, quiet
                                )

                                if deintronize == "on":
                                    log_and_print(
                                        "scanning chiseled sequence for cryptic "
                                        "introns", log_file
                                    )
                                    log_and_print(
                                        "this may take a while...", log_file
                                    )

                                    intron, dico_donor, dico_acceptor = check_intron(
                                        chiseled_sequence
                                    )

                                    if intron:
                                        log_and_print(
                                            f"{current_seq_id} -- possible cryptic "
                                            f"intron found, retrying\n", log_file
                                        )
                                        if dico_donor:
                                            (donor_seqs, donor_locs_left,
                                             donor_locs_right) = intron_seqs_to_avoid(
                                                dico_donor
                                            )
                                        if dico_acceptor:
                                            (acceptor_seqs, acceptor_locs_left,
                                             acceptor_locs_right) = intron_seqs_to_avoid(
                                                dico_acceptor
                                            )

                                        intron_constraints.extend(
                                            donor_seqs + acceptor_seqs
                                        )
                                        loc_constraints_left.extend(
                                            donor_locs_left + acceptor_locs_left
                                        )
                                        loc_constraints_right.extend(
                                            donor_locs_right + acceptor_locs_right
                                        )

                                        original_sequence = chiseled_sequence
                                        last_try += 1
                                        continue

                            else:
                                # Add 5' and 3' ends to the chiseled sequence
                                chiseled_sequence = (
                                    base_5p_end + original_sequence + base_3p_end
                                )

                        except Exception as e:
                            log_and_print(
                                f"{current_seq_id} FAILED chisel.\n"
                                f"Exception occurred: {e}\n", log_file
                            )
                            log_and_print("Retrying...", log_file)
                        last_try = max_attempts

                    # Copy the chisel to "original_sequence" for deintronize mode
                    original_sequence = chiseled_sequence

                    # Add 5' and 3' ends to the chiseled sequence
                    chiseled_sequence = (
                        base_5p_end + original_sequence + base_3p_end
                    )

                if mode == 'no_mods':
                    chiseled_sequence = (
                        base_5p_end + original_sequence + base_3p_end
                    )
                    log_and_print(
                        "\nHinging without chiseling sequence (no sequence changes)",
                        log_file
                    )

                if mode != 'no_hinge':

                    seg_num = math.ceil(
                        len(chiseled_sequence) / (segment_length - radius * 2)
                    )

                    if two_step == "on" and len(chiseled_sequence) > two_step_length:
                        log_and_print("running two-step mode", log_file)
                        two_step_df_all = pd.DataFrame()
                        log_and_print(f"Length: {len(chiseled_sequence)}", log_file, quiet)
                        log_and_print(
                            f"# segments: {seg_num}", log_file, quiet
                        )

                        try:
                            best_solutions_dict = hinge(
                                chiseled_sequence, current_seq_id, two_step_length
                            )

                        except Exception as e:
                            log_and_print(
                                f"{current_seq_id} FAILED hinge.\n"
                                f"Exception occurred: {e}\n", log_file
                            )

                        if best_solutions_dict:
                            for seq_id, solution_data in best_solutions_dict.items():
                                two_step_df_all = solution_data['dataframe']
                                two_step_df_all = pd.concat([two_step_df_all])
                            two_step_df_all['accession'] = accession
                            two_step_df_all['base_id'] = two_step_df_all['Seq_ID']
                            two_step_df_all['Seq_ID'] = (
                                two_step_df_all['Seq_ID'].astype(str) + '_L0_' +
                                two_step_df_all['Fragment_n'].astype(str)
                            )
                            two_step_df_all['Fragment'] = (
                                two_step_5p_end +
                                two_step_df_all['Fragment'].astype(str) +
                                two_step_3p_end
                            )

                            two_step_df = pd.concat([two_step_df, two_step_df_all])

                            best_solutions_dict = {}

                            for index, df in two_step_df_all.iterrows():
                                seq_id = df['Seq_ID']
                                sequence = df['Fragment']
                                seg_num = math.ceil(
                                    len(df['Fragment']) / (segment_length - radius)
                                )
                                log_and_print(
                                    f"Level 1 # segments: {seg_num}", log_file, quiet
                                )

                                best_solutions_dict = hinge(
                                    sequence, seq_id, segment_length
                                )

                            if best_solutions_dict:
                                for seq_id, solution_data in best_solutions_dict.items():
                                    df = solution_data['dataframe']
                                    df['Tries_for_hinges'] = solution_data['Tries_for_hinges']
                                    df['accession'] = accession
                                    df_all = pd.concat([df_all, df])

                    else:
                        log_and_print(f"Length: {len(chiseled_sequence)}", log_file, quiet)
                        log_and_print(
                            f"# segments: {seg_num}", log_file, quiet
                        )

                        try:
                            (best_solution, seq_sets, fidelity_scores,
                             counter_value_for_best_solution) = find_cut_solution(
                                chiseled_sequence, overhang_sets, radius, ext_overhangs,
                                segment_length, n_tries, potapov_data
                            )

                            if best_solution:
                                best_fidelity = max(fidelity_scores.values())
                                log_and_print(
                                    f"\nOH set fidelity: {best_fidelity:.3f}\n", log_file, quiet
                                )

                            if seq_sets:
                                max_fidelity_index = max(
                                    fidelity_scores, key=fidelity_scores.get
                                )
                                solution_data = seq_sets[max_fidelity_index]
                                solutions = solution_data['solutions']
                                counter_value = solution_data['counter']
                                fidelity = fidelity_scores[max_fidelity_index]
                                df = create_fragments_df(
                                    seq_id, solutions, original_sequence,
                                    chiseled_sequence, fidelity, base_5p_end,
                                    base_3p_end, original_cai, chiseled_cai
                                )
                                df['accession'] = accession
                                df['Seq_ID'] = current_seq_id
                                df['Tries_for_hinges'] = counter_value
                                df_all = pd.concat([df_all, df])

                            else:
                                log_and_print(
                                    f"No solutions sets found for {current_seq_id}\n",
                                    log_file
                                )
                        except Exception as e:
                            log_and_print(
                                f"{current_seq_id} FAILED hinge.\n"
                                f"Exception occurred: {e}\n", log_file
                            )

                if mode == 'no_hinge':
                    chiseled_sequence = (
                        base_5p_end + original_sequence + base_3p_end
                    )

                    df = create_chisels_df(
                        accession, current_seq_id, original_sequence, chiseled_sequence,
                        original_cai, chiseled_cai
                    )
                    df_all = pd.concat([df_all, df])
                seed = seed + 1

        yml_file = f"{yml}"

        with open(yml_file, 'r') as f:
            config = yaml.safe_load(f)

        flattened_constraints = {
            f"constraint_{i+1}": item
            for i, item in enumerate(config["constraints"])
        }

        flattened_optimizations = {}
        if "optimizations" in config:
            flattened_optimizations = {
                f"optimization_{i+1}": item
                for i, item in enumerate(config["optimizations"])
            }

        args_dict = vars(updated_defaults)
        args_dict.update(flattened_constraints)
        args_dict.update(flattened_optimizations)

        log_and_print('\n', log_file, quiet)

        with pd.ExcelWriter(f'{ofile}_all_data.xlsx') as writer:
            df_all.to_excel(writer, sheet_name="pre_indexed_data", index=False)

            if two_step == "on" and not two_step_df.empty:
                two_step_df.to_excel(
                    writer, sheet_name="step_1_fragments", index=False
                )

            params_df = pd.DataFrame.from_dict(
                args_dict, orient='index', columns=['Value']
            )
            params_df.to_excel(writer, sheet_name="run_parameters")

        write_fasta(df_all, f'{ofile}_designed_seqs.fasta')

        read_log_and_identify_failures(tag)

        if mode != "no_hinge":
            cmd = [
                'Rscript',
                'scripts/paste_primers_cmd.R',
                '--input_pops', f'{ofile}_all_data.xlsx',
                '--output_file', ofile,
                '--primer_index', str(primer_index),
                '--input_primers', index_primers,
                '--pcr_5p_cut', str(pcr_5p_cut),
                '--pcr_3p_cut', str(pcr_3p_cut),
                '--max_fragments', str(max_fragments),
                '--two_step', 'on' if (two_step == "on" and not two_step_df.empty) else 'off',
                '--run_type', 'cds'
            ]
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as e:
                log_and_print(f"Failed to run R script: {e}", log_file)

        if repeats > 1:
            print('\nCalculating similarity matrix for repeats\n')
            sequences = read_fasta2(f'{ofile}_designed_seqs.fasta')

            for group, seqs in sequences.items():
                unique_seqs = list(set(seqs))

                print(f"\nGroup: {group}")
                log_file.write(f"\nGroup: {group}\n")
                if len(unique_seqs) < 2:
                    print(
                        f"Skipping group {group} because it has less than 2 "
                        f"unique sequences."
                    )
                    log_file.write(
                        f"Skipping group {group} because it has less than 2 "
                        f"unique sequences.\n"
                    )
                    continue

                distance_matrix = calculate_distance_matrix(unique_seqs)
                identity_matrix = calculate_pairwise_identity(distance_matrix)
                print_matrix(identity_matrix, "Identity Matrix:", log_file, decimals=2)

                average_distance = calculate_average_distance(distance_matrix)
                print(f"Average Pairwise Distance (%): {average_distance:.1f}")
                log_file.write(
                    f"Average Pairwise Distance (%): {average_distance:.1f}\n"
                )

        log_and_print('\n', log_file)
        log_and_print('*' * 80, log_file)
        log_and_print('\n', log_file)
        log_and_print("Done", log_file)

        if quiet=="off":
            nerd_alert()

        log_and_print(
            f"\nCompletion time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n",
            log_file, quiet
        )

        log_file.close()

        if mode == "no_hinge" and tweak_n:
            call_tweaker_2(f"out/{tag}/log.txt", tweak_n)

            with open(f'out/{tag}/cai_max_dif_set.fasta', 'w') as fasta_out, \
                    open(f'out/{tag}/cai_max_diff_table.txt', 'w') as table_out:
                for record in SeqIO.parse(f'out/{tag}/max_diff_subsets.fasta', 'fasta'):
                    codon_table = get_codons_table(species)
                    cai = calculate_cai(str(record.seq), codon_table)
                    record.description += f" CAI:{cai:.3f}"
                    SeqIO.write(record, fasta_out, 'fasta')
                    table_out.write(f'{record.id}\t{cai:.3f}\n')

    except Exception as e:
        with open(log_file_path, "a") as log_file:
            log_and_exit(f"An unexpected error occurred: {e}", log_file, exit_code=99)
        log_file.close()
