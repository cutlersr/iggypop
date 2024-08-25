"""
This program is designed to create barcoded fragments of target coding 
sequences that can be amplified in a single PCR reaction and then assembled 
using Golden Gate cloning. Input sequences are optionally optimized using 
DNAchisel. Depending on the parameters set, you can domesticate sequences 
to remove restriction enzymes, perform codon optimization, prevent hairpins, 
and apply various other optimization tools offered by DNAchisel prior to 
generating the fragmented/barcoded sequences. This program is similar to 
'iggypop.py cds' but uses GenBank formatted files as inputs. The sequence 
optimization parameters are set within the GenBank file using DNAchisel's 
rules. GenBank files can be formatted using 'iggypop.py format'; this will 
add the relevant tags for optimization and sequence protection using a YAML 
file to determine optimization patterns.

Once optimized, sequences are appended with defined 5' and 3' ends to create 
overhangs and other elements needed for Golden Gate cloning. The final target 
sequences are fragmented into segments with high-fidelity overhangs for 
efficient reassembly via Golden Gate cloning. The maximum fragment size is 
set by the `--segment_length` flag, with a default value of 200, which is 
appropriate for synthesizing 250 bp oligos (reserving 50 bp for enzyme sites 
and barcode primers). Fragmentation is done with GoldenHinges, which limits 
overhang selection to predefined sets. To ensure high-fidelity assemblies, 
sets of precomputed high-fidelity overhangs are used as constraints. Multiple 
overhang sets are considered (`--n_tries`), with the highest fidelity set 
reported.

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

The script outputs the CAI (codon adaptation index) values for every ORF/CDS 
in the input file(s).

Hinging parameters can be set via command line or more conveniently with a YAML 
file. Command-line arguments override YAML settings. As noted above, sequence 
optimization parameters must be specified in the input GenBank file.

**Command line options**:
- `--i` (str): Input FASTA file path for sequences to optimize and fragment.
- `--o` (str): Output file stem, saved in the `out/` folder.
- `--yml` (str): YAML file for run parameters (default: `yaml/gb_mcu.yml`). 
  Command-line arguments override YAML settings.
- `--mode` (str): Operation mode:
  - `chisel`: Modifies the input sequence (default).
  - `no_mods`: Only hinges and barcodes fragments.
  - `no_hinge`: Chisels sequences without hinging/barcoding.
- `--reports` (bool): Enables DNAchisel reports.

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
"""

from datetime import date
import random
import numpy
import math
import subprocess
import warnings
import pandas as pd
from headers import *
from initialization import *
from chisel_hinge import *
from pop_helpers import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from goldenhinges import (
    OverhangsSelector, list_overhangs, gc_content
)
from dnachisel import (
    DnaOptimizationProblem
)

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
            f"\nOH set fidelity: {best_overall_solution:.3f}\n", log_file
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

    run_type = "gb"
    gb = True
    tag, ofile, log_file_path, updated_defaults = initialize(run_type)
    globals().update(vars(updated_defaults))
    log_file = open(log_file_path, "a")

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

    if tweak_n and mode != "no_hinge":
        print("Tweaking only works in no_hinge mode; run with '--mode no_hinge'.")
        sys.exit()

    print_header_gb(tag)

    # Copy the input files and analysis scripts to the results folder
    gb_file_path = f'{i}'
    yml_file_path = f'{yml}'
    script_paths = [
        "iggypop.py", "iggypop/iggypop_gb.py", "iggypop/pop_helpers.py",
        "iggypop/chisel_hinge.py", "iggypop/initialization.py"
    ]
    output_fasta = f"out/{tag}/{tag}_designed_seqs.fasta"
    output_genbank = f"out/{tag}/{tag}_designed_seqs.gb"

    copy_analysis_assets(gb_file_path, yml_file_path, tag, script_paths)

    with open(f"out/{tag}/log.txt", "a") as log_file:

        log_file.write(f"Log file for: {tag}\n")
        log_file.write(
            f"Start time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        )
        log_file.write(f"Command line: {' '.join(sys.argv)}\n")

        with open(i, "r") as handle:

            df_all = pd.DataFrame()
            two_step_df = pd.DataFrame()

            for record in SeqIO.parse(handle, "genbank"):
                original_sequence = str(record.seq)
                seq_id = record.id  # Extract sequence ID from the record

                for repeat in range(repeats):

                    try:
                        accession = seq_id
                        current_seq_id = f"{seq_id}.{repeat + 1}"
                        best_solutions_dict = {}
                        best_solution = None
                        best_fidelity = 0
                        original_cai = None  
                        chiseled_cai = None  

                        log_and_print(f"\n\nProcessing {accession}\n", log_file)
                        log_and_print('*' * 80, log_file)
        
                        log_and_print(
                            f'>{seq_id}\n{original_sequence}', log_file, quiet
                        )

                        seed = seed if seed else random.randint(0, 2**32 - 1)
                        numpy.random.seed(seed)
                        log_and_print(f'seed: {seed}', log_file, quiet)

                        if mode != 'no_mods':
                            try:

                                problem = DnaOptimizationProblem.from_record(record)

                                log_and_print(
                                    "\nBefore optimization:\n", log_file, quiet
                                )
                                log_and_print(
                                    problem.constraints_text_summary(), log_file, quiet
                                )
                                log_and_print(
                                    problem.objectives_text_summary(), log_file, quiet
                                )

                                problem.resolve_constraints(final_check=True)
                                problem.max_random_iters = 1000

                                results_path = (
                                    f'out/{tag}/reports/{current_seq_id}'
                                )

                                if reports:
                                    problem.optimize_with_report(target=results_path)
                                else:
                                    problem.optimize()

                                log_and_print(
                                    "\nAfter optimization:\n", log_file, quiet
                                )
                                log_and_print(
                                    problem.constraints_text_summary(), log_file, quiet
                                )
                                log_and_print(
                                    problem.objectives_text_summary(), log_file, quiet
                                )

                            except Exception as e:
                                log_and_print(f"An error occurred: {e}", log_file)

                            # Update original_sequence before appending the ends
                            original_sequence = problem.sequence

                            chiseled_sequence = str(
                                base_5p_end + original_sequence + base_3p_end
                            )
                            log_and_print(
                                f'>{current_seq_id}_chisel\n{chiseled_sequence}',
                                log_file, quiet
                            )

                        else:                
                            chiseled_sequence = str(
                                base_5p_end + original_sequence + base_3p_end
                            )
                            log_and_print(
                                f"\nHinging without chiseling sequence "
                                f"(no sequence changes)", log_file
                            )
                            log_and_print(
                                f"Length: {len(chiseled_sequence)}", log_file, quiet
                            )
                            log_and_print(
                                f"# segments: {math.ceil(len(chiseled_sequence) / (segment_length - radius * 2))}",
                                log_file,
                                quiet
                            )

                        if mode != 'no_hinge':

                            seg_num = math.ceil(
                                len(chiseled_sequence) / (segment_length - radius * 2)
                            )

                            if (two_step == "on") and len(chiseled_sequence) > two_step_length:
                                log_and_print("running two-step mode", log_file)
                                two_step_df_all = pd.DataFrame()

                                log_and_print(
                                    f"Length: {len(chiseled_sequence)}", log_file, quiet
                                )
                                log_and_print(f"# segments: {seg_num}", log_file, quiet)
                                
                                try:
                                    best_solutions_dict = hinge(
                                        chiseled_sequence, current_seq_id, two_step_length
                                    )
                                
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
                                        seg_num = math.ceil(len(df['Fragment']) / 
                                                            (segment_length - radius))
                                        log_and_print(
                                            f"Level 1 # segments: {seg_num}", log_file
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

                                except Exception as e:
                                    log_and_print(
                                        f"{current_seq_id} FAILED hinge.\n"
                                        f"Exception occurred: {e}\n", log_file
                                    )


                            else:
                                log_and_print(
                                    f"Length: {len(chiseled_sequence)}", log_file, quiet
                                )
                                log_and_print(f"# segments: {seg_num}", log_file, quiet)

                                try:
                                    (best_solution, seq_sets, fidelity_scores,
                                     counter_value_for_best_solution) = find_cut_solution(
                                        chiseled_sequence, overhang_sets, radius, ext_overhangs,
                                        segment_length, n_tries, potapov_data
                                    )

                                    if best_solution:
                                        best_fidelity = max(fidelity_scores.values())
                                        log_and_print(
                                            f"\nOH set fidelity: {best_fidelity:.3f}\n", 
                                            log_file, quiet
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
                                            base_3p_end, original_cai, chiseled_cai,
                                            gb
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
                            chiseled_sequence = str(
                                base_5p_end + original_sequence + base_3p_end
                            )
                            df = create_chisels_df(
                                accession, current_seq_id, original_sequence, chiseled_sequence,
                                original_cai, chiseled_cai, gb
                            )
                            df['accession'] = accession
                            df_all = pd.concat([df_all, df])

                        # Define the length of base_5p_end
                        base_5p_length = len(base_5p_end)

                        # Update original_sequence before appending the ends
                        original_sequence = problem.sequence if mode != 'no_mods' else str(record.seq)

                        # Append ends to the chiseled sequence
                        chiseled_sequence = str(base_5p_end + original_sequence + base_3p_end)

                        seed += 1

                        with open(output_fasta, 'a') as fasta_file:
                            fasta_file.write(f">{current_seq_id}\n{chiseled_sequence}\n")

                        # Create a new SeqRecord with the chiseled sequences
                        chiseled_seq_record = SeqRecord(
                            Seq(chiseled_sequence), 
                            id=current_seq_id,
                            name=accession,
                            description=record.description,
                            annotations=record.annotations
                        )

                        # Get the current date
                        current_date = datetime.datetime.today().strftime("%d-%b-%Y").upper()

                        # Ensure the sequence length is correct in the LOCUS line
                        chiseled_seq_record.annotations["topology"] = "linear"
                        chiseled_seq_record.annotations["data_file_division"] = "UNA"
                        chiseled_seq_record.annotations["date"] = current_date

                        # Adjust the feature locations
                        adjusted_features = adjust_feature_locations(record.features, base_5p_length)
                        chiseled_seq_record.features = adjusted_features

                        # Update the length annotation directly if needed 
                        # (though length should be inferred correctly)
                        chiseled_seq_record.annotations["sequence_length"] = len(chiseled_seq_record.seq)

                        # Write to the GenBank file
                        output_genbank = f"out/{tag}/{tag}_designed_seqs.gb"
                        with open(output_genbank, 'a') as gb_file:
                            SeqIO.write(chiseled_seq_record, gb_file, "genbank")

                        # Call the report_gb_cai function
                        input_df, output_df = report_gb_cai(gb_file_path, output_genbank, species)

                    except Exception as e:
                        # Print the error message
                        log_and_print(f"An error occurred: {e}", log_file)

            output_path = f'{ofile}_CAIs.txt'
            # Save the CAI summaries to a file
            save_cai_summary_to_file(input_df, output_df, output_path)
            input_df = input_df.drop(columns=['Sequence'])
            output_df = output_df.drop(columns=['Sequence'])
            log_and_print(f'original CAI values for CDSs \n {input_df}\n', log_file)
            log_and_print(f'output CAI values for CDSs \n {output_df}\n', log_file)
                
            # Convert command-line arguments to a dictionary
            dict = vars(updated_defaults)

            # Create output Excel file and then add primers and format output with R
#            if mode != 'no_hinge':
#                # Drop unnecessary columns for gb run
#                df_all = df_all.drop(columns=['Original_CAI', 'Chiseled_CAI'])

            # Create an ExcelWriter object
            with pd.ExcelWriter(f'{ofile}_all_data.xlsx') as writer:
                # Write the main dataframe to the first sheet
                df_all.to_excel(writer, sheet_name="pre_indexed_data", index=False)
                if two_step and not two_step_df.empty:
                    two_step_df.to_excel(writer, sheet_name="step_1_fragments", index=False)
                # Convert the merged dict to a DataFrame and write to a new sheet
                params_df = pd.DataFrame.from_dict(dict, orient='index', columns=['Value'])
                params_df.to_excel(writer, sheet_name="run_parameters")

            if repeats > 1:
                try:
                    log_and_print('\nCalculating similarity matrix for repeats\n', log_file)
                    sequences = read_fasta2(f'{ofile}_designed_seqs.fasta')

                    for group, seqs in sequences.items():
                        # Filter unique sequences
                        unique_seqs = list(set(seqs))

                        log_and_print(f"\nGroup: {group}", log_file)

                        if len(unique_seqs) < 2:
                            log_and_print(
                                f"Skipping group {group} because it has less than 2 "
                                f"unique sequences.", log_file
                            )
                            continue

                        distance_matrix = calculate_distance_matrix(unique_seqs)
                        identity_matrix = calculate_pairwise_identity(distance_matrix)
                        print_matrix(identity_matrix, "Identity Matrix:", log_file, decimals=2)

                        average_distance = calculate_average_distance(distance_matrix)
                        log_and_print(
                            f"Average Pairwise Distance (%): {average_distance:.1f}", 
                            log_file
                        )

                except Exception as e:
                    # Print the error message
                    log_and_print(f"An error occurred: {e}", log_file)

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
                    '--run_type', str(run_type)
                ]
                # Execute the R script with error handling
                try:
                    subprocess.run(cmd, check=True)
                except subprocess.CalledProcessError as e:
                    log_and_print(f"Failed to run R script: {e}", log_file)

            log_file.write(
                f"\nCompletion time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
            )

            log_file.close()

            if mode == "no_hinge" and tweak_n:
                call_tweaker_2(f"out/{tag}/log.txt", tweak_n)

            read_log_and_identify_failures(tag)

            if quiet=="off":
                nerd_alert()
