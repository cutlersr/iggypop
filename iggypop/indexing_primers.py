import os
import tempfile
import random
import re
from multiprocessing import Pool, cpu_count
import multiprocessing
import pandas as pd
import argparse
import primer3
import glob
from pathlib import Path
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio.Seq import Seq
from concurrent.futures import ProcessPoolExecutor
import datetime

"""
This script is used to design primer pairs that act as indexing barcodes 
to allow amplification of specific targets from oligo pools. It generates 
sets of primers of a specified length and Tm using Primer3, and then 
filters them using MFEPrimer to select the primer pairs (i) lacking 
hairpins (ii) not prone to primer dimers and (iii) with the lowest 
likelihood of cross-hybridization between primers (specified by 
--first_cut percentile for the calculated off-target binding scores.)

Experimental validation of ~400 18 bp primer pairs (60 oC Tm) designed 
with this pipeline show ~99% successful amplification of targets from 
complex oligo pools. A set of 10K primers made with this pipeline is 
provided in the data folder.

Primer Design:
   - Generates random DNA sequences of a specified length and GC content 
     lacking specified restriction sites & sequences; these are then fed 
     to Primer3 to select optimal PCR primer pairs for the sequence with a 
     specified primer length and melting temperature (Tm).

Primer Filtering:
   - Analyzes the designed primers using MFEPrimer to identify potential 
     issues such as cross-hybridization, hairpin formation, and 
     dimerization. Filters the primer pairs to identify and exclude 
     primers that could bind to E. coli and other common contaminants 
     DNAs using MFEPriner and a database of contaminant DNAs; primers 
     are ranked based on their binding scores for each primer; pairs 
     that generate potential amplicons are removed.

Parameters:
- `num_sequences`: Number of DNA sequences to generate.
- `sequence_length`: Length of each DNA sequence.
- `gc_content`: Desired GC content of the sequences.
- `opt_size`: Optimal size of the primers.
- `opt_tm`: Optimal melting temperature of the primers.
- `max_size`: Maximum size of the primers.
- `min_size`: Minimum size of the primers.
- `cores`: Number of CPU cores to use for parallel processing.
- `no_coli`: Flag to exclude search against E. coli contaminants.

Outputs to `out/primers` directory.
MFEprimer executable needs to be in the data directory; a linux version 
is included in the repo.
"""

# Sequences to exclude from primers:
# BsaI, BbsI, BfuAO,BsmBI, BspMI, Esp3I, PaqCI, SapI, EcoRI, HindIII
# BamHI, PstI, SalI, KpnI, XbaI, NotI, SacI, XhoI, SpeI, NcoI
forbidden_sequences = [
    "GGTCTC", "GAGACC", "GAAGAC", "GTCTTC", "ACCTGC", "GCAGGT",
    "CGTCTC", "GAGACG", "ACCTGC", "GCAGGT", "CGTCGC", "GCGACG",
    "CACCTGC", "GCAGGTG", "GCTCTCTTC", "GAAGAGAGC", "GAATTC",
    "AAGCTT", "GGATCC", "CTGCAG", "GTCGAC", "GGTACC", "TCTAGA",
    "GCGGCCGC", "GAGCTC", "CTCGAG", "ACTAGT", "CCATGG"
]


def generate_dna_sequence(sequence_length, gc_content):
    num_gc = round(sequence_length * gc_content)
    num_at = sequence_length - num_gc

    num_g = num_c = num_gc // 2
    num_a = num_t = num_at // 2

    bases = ['G'] * num_g + ['C'] * num_c + ['A'] * num_a + ['T'] * num_t

    for _ in range(10000):
        random.shuffle(bases)
        sequence = ''.join(bases)
        if not any(fs in sequence for fs in forbidden_sequences):
            return sequence
    return None


def worker(args):
    return generate_dna_sequence(*args)


def generate_dna_sequences(num_sequences, sequence_length, gc_content, cores):
    pool = Pool(cores)
    sequences = pool.map(
        worker,
        [(sequence_length, gc_content) for _ in range(num_sequences)]
    )
    pool.close()
    pool.join()

    sequences_df = pd.DataFrame({
        'sequence': sequences,
        'length': sequence_length,
        'gc_content': gc_content
    })

    return sequences_df


def reverse_complement(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complemented = (complement[nuc] for nuc in dna_sequence.upper())
    return ''.join(reversed(list(complemented)))


def generate_pcr_primers(sequence, opt_size=18, opt_tm=60.0,
                         max_size=18, min_size=18):
    primers = primer3.bindings.design_primers(
        {
            'SEQUENCE_ID': 'seq1',
            'SEQUENCE_TEMPLATE': sequence,
        },
        {
            'PRIMER_OPT_SIZE': opt_size,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_MIN_SIZE': min_size,
            'PRIMER_MAX_SIZE': max_size,
            'PRIMER_OPT_TM': opt_tm,
            'PRIMER_OPT_GC_PERCENT': 55.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 40.0,
            'PRIMER_MAX_GC': 60.0,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_SELF_ANY': 8,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 8,
            'PRIMER_PAIR_MAX_COMPL_END': 8.0,
            'PRIMER_MAX_END_GC': 3,
            'PRIMER_MAX_END_STABILITY': 7,
        }
    )
    return primers


def main(num_sequences, sequence_length, gc_content, opt_size, opt_tm,
         max_size, min_size, mfeprimer_path='./data/mfeprimer-3.3.1-linux-386',
         cores=11):
    sequences_df = generate_dna_sequences(
        num_sequences, sequence_length, gc_content, cores
    )
    primer_data = []

    print(
        f'Generating {opt_size} bp indexing primers with ~Tm of {opt_tm} '
        f'using Primer3'
    )

    for index, sequence in enumerate(sequences_df['sequence']):
        primers = generate_pcr_primers(
            sequence, opt_size, opt_tm, max_size, min_size
        )
        if 'PRIMER_LEFT_0_SEQUENCE' in primers and \
                'PRIMER_RIGHT_0_SEQUENCE' in primers:
            primer_data.append({
                'name': f'sequence_{index}',
                'left_seq': primers['PRIMER_LEFT_0_SEQUENCE'],
                'right_seq': primers['PRIMER_RIGHT_0_SEQUENCE'],
                'left_tm': primers['PRIMER_LEFT_0_TM'],
                'right_tm': primers['PRIMER_RIGHT_0_TM'],
                'left_gc_pct': primers['PRIMER_LEFT_0_GC_PERCENT'],
                'right_gc_pct': primers['PRIMER_RIGHT_0_GC_PERCENT'],
                'left_self_any_th': primers['PRIMER_LEFT_0_SELF_ANY_TH'],
                'right_self_any_th': primers['PRIMER_RIGHT_0_SELF_ANY_TH'],
                'left_self_end_th': primers['PRIMER_LEFT_0_SELF_END_TH'],
                'right_self_end_th': primers['PRIMER_RIGHT_0_SELF_END_TH'],
                'left_hairpin_th': primers['PRIMER_LEFT_0_HAIRPIN_TH'],
                'right_hairpin_th': primers['PRIMER_RIGHT_0_HAIRPIN_TH'],
                'left_end_stability': primers['PRIMER_LEFT_0_END_STABILITY'],
                'right_end_stability': primers['PRIMER_RIGHT_0_END_STABILITY'],
            })
        else:
            print(f"No primers found for sequence_{index}")

    primers_df = pd.DataFrame(primer_data)

    primers_df = primers_df[
        (primers_df['left_seq'].apply(len) == opt_size) &
        (primers_df['right_seq'].apply(len) == opt_size)
    ]

    print(
        'Scanning indexing primers for possible cross-hybridization and other '
        'problems using mfeprimer'
    )

    create_target_seqs_for_mfe(primers_df)
    subprocess.run(
        [mfeprimer_path, "index", "-i", "out/primers/all_pri.fa"],
        check=True
    )

    return primers_df


def create_target_seqs_for_mfe(primers_df, stuffer_length=150,
                               mfeprimer_path='./data/mfeprimer-3.3.1-linux-386'):
    all_primers_lines = []

    for index, row in primers_df.iterrows():
        right_seq_revcomp = str(
            Seq(row['right_seq']).reverse_complement()
        )

        combination = f"{row['left_seq']}{'A' * stuffer_length}{right_seq_revcomp}"

        all_primers_lines.append({
            'name': f">{row['name']}",
            'sequence': combination
        })

    new_primers_df = pd.DataFrame(all_primers_lines)

    output_dir = "out/primers"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_file = os.path.join(output_dir, "all_pri.fa")

    with open(output_file, 'w') as fasta_file:
        for index, row in new_primers_df.iterrows():
            fasta_file.write(f"{row['name']}\n{row['sequence']}\n")

    return new_primers_df


def parse_mfeprimer_output(output_file):
    with open(output_file, 'r') as file:
        lines = file.readlines()

    no_hairpins = any("No Hairpins found." in line for line in lines)
    no_dimers = any("No dimers found." in line for line in lines)

    amplicon_line = next(
        (line for line in lines if "Descriptions of [" in line), None
    )
    if amplicon_line:
        amplicon_match = re.search(r"\[\s*(\d+)\s*\]", amplicon_line)
        no_amplicons = int(amplicon_match.group(1)) if amplicon_match else pd.NA
    else:
        no_amplicons = pd.NA

    primer_pattern = re.compile(
        r'(set\d+_(?:left|right))\s+([\w-]+)\s+(\d+)\s+(\d+\.\d+)\s+'
        r'(\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+)\s+(\d+)'
    )
    data = []

    for line in lines:
        match = primer_pattern.search(line)
        if match:
            primer_id, sequence, length, gc, tm, dg, binding_plus, \
                binding_minus = match.groups()
            data.append({
                'primer_id': primer_id,
                'sequence': sequence,
                'bp': int(length),
                'gc_pct': float(gc),
                'Tm': float(tm),
                'Dg': float(dg),
                'binding_plus': int(binding_plus),
                'binding_minus': int(binding_minus),
                'hairpins': 0 if no_hairpins else 1,
                'dimers': 0 if no_dimers else 1,
                'amplicons': no_amplicons
            })

    if not data:
        print("No data extracted, check the file and regular expression.")
        return pd.DataFrame()

    df = pd.DataFrame(data)

    if len(df) >= 2:
        binding_numbers_df = pd.DataFrame({
            'F_bind_plus': [df.iloc[0]['binding_plus']],
            'F_bind_minus': [df.iloc[0]['binding_minus']],
            'R_bind_plus': [df.iloc[1]['binding_plus']],
            'R_bind_minus': [df.iloc[1]['binding_minus']],
            'F_Tm': [df.iloc[0]['Tm']],
            'F_Dg': [df.iloc[0]['Dg']],
            'R_Tm': [df.iloc[1]['Tm']],
            'R_Dg': [df.iloc[1]['Dg']],
            'hairpins': [0 if no_hairpins else 1],
            'dimers': [0 if no_dimers else 1],
            'amplicons': [no_amplicons]
        })
        return binding_numbers_df
    else:
        print("Not enough primer data to create binding numbers DataFrame.")
        return pd.DataFrame()


def mfeprimer(seq_num, all_primers, mfeprimer_path):
    temp_file_path = "out/primers/test.fasta"
    out_file_path = "out/primers/out.txt"

    primer_lines = [
        f">set{seq_num}_left",
        all_primers.loc[seq_num, 'left_seq'],
        f">set{seq_num}_right",
        all_primers.loc[seq_num, 'right_seq']
    ]

    with open(temp_file_path, 'w') as file:
        file.write("\n".join(primer_lines))

    subprocess.run(
        [mfeprimer_path, "-i", temp_file_path, "-d", "out/primers/all_pri.fa",
         "-o", out_file_path], check=True
    )

    analysis = parse_mfeprimer_output(out_file_path)
    analysis['name'] = all_primers.loc[seq_num, 'name']

    return analysis


def run_mfeprimer(all_primers,
                  mfeprimer_path='./data/mfeprimer-3.3.1-linux-386'):
    analysis_results = []

    for index in range(len(all_primers)):
        result = mfeprimer(index, all_primers, mfeprimer_path)
        analysis_results.append(result)

    analysis_df = pd.concat(analysis_results, ignore_index=True)

    os.remove("out/primers/test.fasta")
    os.remove("out/primers/out.txt")
    for file in glob.glob("out/primers/all_pri.fa*"):
        os.remove(file)
    return analysis_df


def mfe_coli(seq_num, primers, mfeprimer_path, database):
    primer_lines = [
        f">set{seq_num}_left\n{primers.loc[seq_num, 'left_seq']}",
        f">set{seq_num}_right\n{primers.loc[seq_num, 'right_seq']}"
    ]

    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
        temp_file_path = temp_file.name
        temp_file.write("\n".join(primer_lines))
        temp_file.flush()

    output_file_path = temp_file_path + "_output"
    subprocess.run(
        [mfeprimer_path, "-i", temp_file_path, "-d", database,
         "-o", output_file_path], check=True
    )

    analysis = parse_mfeprimer_output(output_file_path)
    analysis['name'] = primers.loc[seq_num, 'name']

    os.unlink(temp_file_path)
    os.unlink(output_file_path)

    return analysis


def run_mfe_coli_all(primers, database, mfeprimer_path, num_cores=10):
    with ThreadPoolExecutor(max_workers=num_cores) as executor:
        futures = [
            executor.submit(
                mfe_coli, i, primers, mfeprimer_path, database
            ) for i in range(len(primers))
        ]
        results = [future.result() for future in as_completed(futures)]

    all_analysis = pd.concat(results, ignore_index=True)
    return all_analysis


def mfeprimer(seq_num, all_primers, mfeprimer_path, temp_file_suffix):
    temp_file_path = f"out/primers/test_{temp_file_suffix}.fasta"
    out_file_path = f"out/primers/out_{temp_file_suffix}.txt"

    primer_lines = [
        f">set{seq_num}_left",
        all_primers.loc[seq_num, 'left_seq'],
        f">set{seq_num}_right",
        all_primers.loc[seq_num, 'right_seq']
    ]

    with open(temp_file_path, 'w') as file:
        file.write("\n".join(primer_lines))

    subprocess.run(
        [mfeprimer_path, "-i", temp_file_path, "-d", "out/primers/all_pri.fa",
         "-o", out_file_path], check=True
    )

    analysis = parse_mfeprimer_output(out_file_path)
    analysis['name'] = all_primers.loc[seq_num, 'name']

    os.remove(temp_file_path)
    os.remove(out_file_path)

    return analysis


def run_mfeprimer_multicore(all_primers,
                            mfeprimer_path='./data/mfeprimer-3.3.1-linux-386',
                            num_cores=10):
    analysis_results = []

    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        futures = {
            executor.submit(
                mfeprimer, index, all_primers, mfeprimer_path, index
            ): index for index in range(len(all_primers))
        }

        for future in as_completed(futures):
            try:
                result = future.result()
                analysis_results.append(result)
            except Exception as e:
                print(f"Task {futures[future]} generated an exception: {e}")

    if analysis_results:
        analysis_df = pd.concat(analysis_results, ignore_index=True)
    else:
        analysis_df = pd.DataFrame()

    for file in glob.glob("out/primers/all_pri.fa*"):
        os.remove(file)

    return analysis_df


def get_unique_filename(base_path, filename):
    base_name, ext = os.path.splitext(filename)
    full_path = os.path.join(base_path, filename)
    counter = 1

    while os.path.exists(full_path):
        full_path = os.path.join(base_path, f"{base_name}.{counter}{ext}")
        counter += 1

    return full_path


def check_mfeprimer():
    data_dir = './data'
    mfeprimer_glob = os.path.join(data_dir, 'mfeprimer*')
    mfeprimer_versions = glob.glob(mfeprimer_glob)

    if len(mfeprimer_versions) == 0:
        print("MFEprimer is not installed in the data folder.")
        print("Please download MFEprimer using the following commands:")
        print("""
        wget https://github.com/quwubin/MFEprimer-3.0/releases/download/v3.3.1/mfeprimer-3.3.1-linux-386.gz
        gunzip mfeprimer-3.3.1-linux-386.gz
        mv mfeprimer-3.3.1-linux-386 ./data/
        chmod +x ./data/mfeprimer-3.3.1-linux-386
        """)
        exit(1)
    elif len(mfeprimer_versions) == 1 and \
            mfeprimer_versions[0].endswith('mfeprimer-3.3.1-linux-386'):
        print("MFEprimer installed: mfeprimer-3.3.1-linux-386")
    else:
        for mfe_version in mfeprimer_versions:
            if 'mfeprimer-3.3.1-linux-386' in mfe_version:
                print("MFEprimer installed: mfeprimer-3.3.1-linux-386")
            else:
                print(f"Warning: Detected MFEprimer version '{mfe_version}'.")
                print(
                    "Iggypop was not specifically designed with this version, "
                    "but it should work fine."
                )

    return mfeprimer_versions[0]  # Return the detected MFEprimer path


def rescore_existing_file(file_path, mfeprimer_path, num_cores):
    primers_df = pd.read_csv(file_path)

    required_columns = [
        'name', 'F_seq', 'R_seq', 'R_seq_rc', 'F_Tm', 'R_Tm', 'Tm_avg',
        'F_Dg', 'R_Dg', 'Dg_avg', 'pcr_check', 'relocated from'
    ]
    if not all(col in primers_df.columns for col in required_columns):
        raise ValueError(f"The CSV file must contain the following columns: {required_columns}")

    primers_df = primers_df.rename(columns={
        'F_seq': 'left_seq',
        'R_seq': 'right_seq'
    })

    # Create target sequences with stuffer sequences included
    print("Creating target sequences with stuffer sequences included...")

    primers_with_stuffer_df = primers_df.copy()
    primers_with_stuffer_df['sequence'] = primers_with_stuffer_df.apply(
        lambda row: f"{row['left_seq']}{'A' * 150}{row['R_seq_rc']}",
        axis=1
    )

    output_dir = "out/primers"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    all_pri_file = os.path.join(output_dir, "all_pri.fa")
    with open(all_pri_file, 'w') as fasta_file:
        for index, row in primers_with_stuffer_df.iterrows():
            fasta_file.write(f">{row['name']}\n{row['sequence']}\n")

    print(f"FASTA file with stuffer sequences created at {all_pri_file}")

    subprocess.run(
        [mfeprimer_path, "index", "-i", "out/primers/all_pri.fa"],
        check=True
    )

    try:
        primers_analysis_df = run_mfeprimer_multicore(
            primers_with_stuffer_df, mfeprimer_path, num_cores
        )
    except Exception as e:
        print(f"Error during MFEPrimer self-binding analysis: {e}")
        raise

    # Run MFE analysis against the contaminant database
    try:
        primers_coli_df = run_mfe_coli_all(
            primers_with_stuffer_df, all_pri_file, mfeprimer_path, num_cores
        )
    except Exception as e:
        print(f"Error during MFEPrimer contaminant-binding analysis: {e}")
        raise

    try:
        primers_df = primers_df.merge(primers_analysis_df, on='name', how='left')
        primers_df = primers_df.merge(primers_coli_df, on='name', how='left')
    except Exception as e:
        print(f"Error during merging of analysis results: {e}")
        raise

    primers_df.columns = primers_df.columns.str.replace(r'_x$', '', regex=True)
    primers_df.columns = primers_df.columns.str.replace(r'_y$', '', regex=True)

    primers_df = primers_df.loc[:, ~primers_df.columns.duplicated()]

    # Calculate the off-target binding score
    try:
        primers_df['off_target_binding'] = (
            primers_df[['F_bind_plus', 'F_bind_minus', 'R_bind_plus', 'R_bind_minus']].sum(axis=1) - 2
        ).fillna(0)
    except KeyError as e:
        print(f"Error: Column {e} not found in DataFrame. Available columns: {primers_df.columns}")
        raise

    if 'amplicons' in primers_df.columns:
        primers_df = primers_df.drop(columns=['amplicons'])

    # Calculate contaminant_binding and place it at the end
    try:
        primers_df['contaminant_binding'] = (
            primers_df[['F_bind_plus', 'F_bind_minus', 'R_bind_plus', 'R_bind_minus']].sum(axis=1)
        ).fillna(0)
    except KeyError as e:
        print(f"Error: Column {e} not found in DataFrame. Available columns: {primers_df.columns}")
        raise

    primers_df = primers_df.rename(columns={
        'left_seq': 'F_seq',
        'right_seq': 'R_seq'
    })

    # Select and reorder the relevant columns for the final output
    output_columns = [
        'name', 'F_seq', 'R_seq', 'R_seq_rc', 'F_Tm', 'R_Tm', 'Tm_avg',
        'F_Dg', 'R_Dg', 'Dg_avg', 'hairpins', 'dimers', 'off_target_binding',
        'pcr_check', 'relocated from', 'contaminant_binding'
    ]

    for col in output_columns:
        if col not in primers_df.columns:
            raise KeyError(f"Column '{col}' not found in DataFrame. Available columns: {primers_df.columns}")

    primers_df = primers_df[output_columns]

    output_path = get_unique_filename(
        os.path.dirname(file_path), "rescored_" + os.path.basename(file_path)
    )
    primers_df.to_csv(output_path, index=False)
    print(f"Rescored data saved to {output_path}")
    return primers_df

if __name__ == '__main__':
    mfeprimer_path = check_mfeprimer()

    parser = argparse.ArgumentParser(
        description='Generate DNA sequences and their PCR primers or rescore an existing CSV file.'
    )
    parser.add_argument(
        '--num_sequences', type=int, default=10,
        help='Number of sequences to generate'
    )
    parser.add_argument(
        '--sequence_length', type=int, default=800,
        help='Length of each sequence'
    )
    parser.add_argument(
        '--gc_content', type=float, default=0.5,
        help='GC content ratio'
    )
    parser.add_argument(
        '--opt_size', type=int, default=18,
        help='Optimal size of the primer.'
    )
    parser.add_argument(
        '--opt_tm', type=float, default=60.0,
        help='Optimal melting temperature of the primer.'
    )
    parser.add_argument(
        '--max_size', type=int, default=18,
        help='Maximum size of the primer.'
    )
    parser.add_argument(
        '--min_size', type=int, default=18,
        help='Minimum size of the primer.'
    )
    parser.add_argument(
        '--cores', type=int, default=11,
        help='Number of cores to use.'
    )
    parser.add_argument(
        '--no_coli', action='store_true', default=False,
        help='Exclude default search for amplicons of coli and other '
             'contaminants'
    )
    parser.add_argument(
        '--first_cut', type=float, default=0.3,
        help='precentile of primers to retain at first filtering; default '
             '≤30th percentile; lowest binding scores.'
    )
    parser.add_argument(
        '--second_cut', type=float, default=0.5,
        help='percentile of primers to retain on second contaminant filtering; '
             'default ≤50th binding scores.'
    )
    parser.add_argument(
        '--rescore', type=str,
        help='Path to an existing CSV file to rescore'
    )

    args = parser.parse_args()

    if args.rescore:
        rescore_existing_file(args.rescore, mfeprimer_path, args.cores)

    else:
        print(f"Start time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

        primers_df = main(
            args.num_sequences, args.sequence_length, args.gc_content,
            args.opt_size, args.opt_tm, args.max_size, args.min_size,
            mfeprimer_path=mfeprimer_path, cores=args.cores
        )

        # Run MFEPrimer analysis
        primer_analysis_df = run_mfeprimer_multicore(
            primers_df, mfeprimer_path, args.cores
        )
        primers_merged = pd.merge(primers_df, primer_analysis_df, on="name")
        primers_merged['off_target_binding'] = (
            primers_merged['F_bind_plus'] + primers_merged['F_bind_minus'] +
            primers_merged['R_bind_plus'] + primers_merged['R_bind_minus']
        )

        # Filter out the bottom nth quantile of records based on 
        # 'off_target_binding' & eliminate primers with hairpins or propensity to 
        # form dimers
        cutoff = primers_merged['off_target_binding'].quantile(args.first_cut)
        primers_merged = primers_merged[
            primers_merged['off_target_binding'] <= cutoff
        ]
        primers_merged.sort_values(by=['off_target_binding'], inplace=True)
        primers_merged['off_target_binding'] = primers_merged[
                                                   'off_target_binding'] - 2
        primers_merged[
            (primers_merged['amplicons'] == 1) &
            (primers_merged['hairpins'] == 0) &
            (primers_merged['dimers']) == 0
        ]

        output_dir = "out/primers"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        all_primers_output_path = get_unique_filename(
            output_dir, "all_primers_post_mfe.csv"
        )
        primers_merged.to_csv(all_primers_output_path, index=False)

        primers_df = pd.read_csv(all_primers_output_path)
        database = "data/coli_other_contaminants.fasta"

        if not args.no_coli:
            print(
                'Removing primers with cross-hybridization with E. coli and other '
                'contaminating DNAs'
            )
            all_analysis = []
            all_analysis = run_mfe_coli_all(
                primers_df, database, mfeprimer_path, args.cores
            )

            clean_primers = all_analysis[
                (all_analysis['amplicons'] == 0) &
                (all_analysis['hairpins'] == 0) &
                (all_analysis['dimers'] == 0)
            ].copy()

            clean_primers['contaminant_binding'] = (
                clean_primers['F_bind_plus'] + clean_primers['F_bind_minus'] +
                clean_primers['R_bind_plus'] + clean_primers['R_bind_minus']
            )
            clean_primers['Tm_avg'] = (
                clean_primers['F_Tm'] + clean_primers['R_Tm']
            ) / 2
            clean_primers['Dg_avg'] = (
                clean_primers['F_Dg'] + clean_primers['R_Dg']
            ) / 2

            columns_to_read = [
                'name', 'hairpins', 'dimers', 'amplicons', 'off_target_binding'
            ]
            off_target_df = pd.read_csv(
                all_primers_output_path, usecols=columns_to_read
            )

            clean_primers = clean_primers.merge(off_target_df, on='name', how='left')

            clean_primers.drop(
                columns=['hairpins_x', 'dimers_x', 'amplicons_x'], inplace=True
            )
            clean_primers.rename(
                columns={
                    'hairpins_y': 'hairpins', 'dimers_y': 'dimers',
                    'amplicons_y': 'amplicons'
                }, inplace=True
            )

            clean_primers.sort_values(
                by=['off_target_binding', 'Tm_avg'], inplace=True
            )
            clean_primers = clean_primers.merge(
                primers_df[['name', 'left_seq', 'right_seq']], on='name'
            )

            required_columns = [
                'name', 'left_seq', 'right_seq', 'F_Tm', 'R_Tm', 'Tm_avg', 'F_Dg',
                'R_Dg', 'Dg_avg', 'contaminant_binding', 'hairpins', 'dimers',
                'amplicons', 'off_target_binding'
            ]
            missing_columns = [
                col for col in required_columns if col not in clean_primers.columns
            ]
            if missing_columns:
                print(f"Missing columns after merge: {missing_columns}")
            else:
                clean_primers = clean_primers[required_columns]

            clean_primers.rename(
                columns={'left_seq': 'F_seq', 'right_seq': 'R_seq'}, inplace=True
            )
            clean_primers = clean_primers.drop_duplicates(
                subset=['F_seq']
            ).drop_duplicates(subset=['R_seq'])
            clean_primers.reset_index(drop=True, inplace=True)
            clean_primers.index += 1

            cutoff = clean_primers['contaminant_binding'].quantile(args.second_cut)
            clean = clean_primers[
                clean_primers['contaminant_binding'] <= cutoff
            ]

            clean = clean.sort_values(
                by=['off_target_binding', 'Tm_avg'],
                ascending=[True, False]
            )
            clean['R_seq_rc'] = clean['R_seq'].apply(reverse_complement)

            clean['name'] = ['set' + str(i) for i in range(1, len(clean) + 1)]

            col_order = [
                'name', 'F_seq', 'R_seq', 'R_seq_rc', 'F_Tm', 'R_Tm', 'Tm_avg',
                'F_Dg', 'R_Dg', 'Dg_avg', 'contaminant_binding', 'hairpins',
                'dimers', 'amplicons', 'off_target_binding'
            ]
            clean = clean[col_order]

            final_primers_output_path = get_unique_filename(
                output_dir, "final_primers.csv"
            )
            clean.to_csv(final_primers_output_path, index=False)

        print(f'Done')
        print(f'Your primers are located in the out/primers folder')
        print(f"Stop time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
