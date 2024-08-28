import os
import sys
import re
import shutil
import time
import datetime
import random
import subprocess
from typing import Tuple, List
import numpy as np
import textwrap
import pandas as pd
import math
import argparse
import yaml
from headers import *
import copy
import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import BiopythonParserWarning
from copy import deepcopy
from goldenhinges import (OverhangsSelector)
from dnachisel import (
    DnaOptimizationProblem, 
    CodonOptimize, 
    AvoidPattern, 
    AvoidHairpins, 
    EnforceGCContent, 
    EnforceTranslation, 
    EnforceChanges, 
    UniquifyAllKmers, 
    AvoidChanges,
    AvoidMatches,
    AvoidBlastMatches,
    AvoidStopCodons,
    EnforceChanges,
    EnforcePatternOccurence)
from python_codon_tables import get_codons_table
from Bio.Align import PairwiseAligner

# Ignore warnings
warnings.simplefilter(action='ignore', category=Warning)


def create_fragments_df(
    seq_id, solutions, original_sequence, chiseled_sequence, fidelity,
    base_5p_end, base_3p_end, original_cai, chiseled_cai, gb=False
):
    """
    Create a dataframe containing fragment information.

    Parameters:
    seq_id (str): Identifier for the sequence.
    solutions (list): List of solutions for the sequence.
    original_sequence (str): Original DNA sequence.
    chiseled_sequence (str): Chiseled DNA sequence.
    fidelity (float): Fidelity score for the sequence.
    base_5p_end (str): 5' end of the sequence.
    base_3p_end (str): 3' end of the sequence.
    original_cai (float): CAI for the original sequence.
    chiseled_cai (float): CAI for the chiseled sequence.
    gb (bool): Flag indicating whether to include protein translations.

    Returns:
    pd.DataFrame: Dataframe containing the fragment information.
    """
    df = report_fragments(seq_id, solutions, chiseled_sequence)
    df['Seq_ID'] = seq_id
    df['Internal_overhangs'] = ", ".join([o["sequence"] for o in solutions])
    df['Fidelity'] = fidelity
    df['Original_sequence'] = original_sequence
    df['Chiseled_sequence'] = chiseled_sequence
    df['Original_CAI'] = original_cai
    df['Chiseled_CAI'] = chiseled_cai
    df['Changes'] = df.apply(
        lambda row: highlight_changes(
            row['Original_sequence'], row['Chiseled_sequence'], base_5p_end
        ), axis=1
    )
    if not gb:
        df['Original_protein'] = str(Seq(original_sequence).translate())
        df['Chiseled_protein'] = str(
            Seq(chiseled_sequence[
                len(base_5p_end):(len(chiseled_sequence) - len(base_3p_end))
            ]).translate()
        )
    return df

def report_fragments(seq_id, solution, sequence):
    """
    Report the fragments of a sequence after hinging it.

    Parameters:
    seq_id (str): Identifier for the sequence.
    solution (list): List of solutions for the sequence.
    sequence (str): The DNA sequence.

    Returns:
    pd.DataFrame: Dataframe containing the fragment information.
    """
    L = len(sequence)
    overhang_length = len(solution[0]["sequence"])
    if solution[0]["location"] != 0:
        solution = [{"location": 0, "sequence": sequence[:overhang_length]}] + solution
    if solution[-1]["location"] != L - overhang_length:
        solution = solution + [
            {
                "location": L - overhang_length,
                "sequence": sequence[L - overhang_length:],
            }
        ]
    data = []
    for i, (o1, o2) in enumerate(zip(solution, solution[1:])):
        start, end = o1["location"], o2["location"] + len(o2["sequence"])
        fragment = sequence[start:end]
        data.append([
            seq_id, i+1, fragment, len(fragment), o1['sequence'], o2['sequence']
        ])
    df = pd.DataFrame(
        data, columns=[
            'Seq_ID', 'Fragment_n', 'Fragment', 'Fragment_length',
            'Left_Overhang', 'Right_Overhang'
        ]
    )

    return df


def highlight_changes(original, chiseled, base_5p_end):
    """
    Highlight the changes between the original and chiseled sequences
    with lower case to highlight changes.


    Parameters:
    original (str): Original DNA sequence.
    chiseled (str): Chiseled DNA sequence.
    base_5p_end (str): 5' end of the sequence.

    Returns:
    str: String highlighting the changes.
    """
    offset = len(base_5p_end)
    original_offset = " " * offset + original
    chiseled_offset = chiseled

    changes = ""
    for o, c in zip(original_offset, chiseled_offset):
        if o == c:
            changes += o.upper()
        else:
            changes += c.lower()
    return changes


def check_orf(
    original_sequence, seq_id, allowed_chars, segment_length, log_file=None
):
    """
    Check if a sequence is an open reading frame (ORF).

    'iggypop.py cds' mode only works with valid ORFs

    Parameters:
    original_sequence (str): The DNA sequence to check.
    seq_id (str): Identifier for the sequence.
    allowed_chars (set): Set of allowed characters in the sequence.
    segment_length (int): Minimum length of the sequence segment.
    log_file (file): Optional log file for logging messages.

    Returns:
    bool: True if the sequence is a complete ORF, False otherwise.
    """
    if (
        len(original_sequence) % 3 != 0 or 
        not set(original_sequence).issubset(allowed_chars) or 
        len(original_sequence) < segment_length or 
        original_sequence[:3] not in ["ATG", "GTG", "TTG", "CTG"]
    ):
        log_message = (
            f"\nSkipping sequence {seq_id} because it's not a complete ORF,\n"
            "doesn't start with ATG/GTG/TTG/CTG, contains invalid characters,\n"
            "or is not long enough.\n"
        )
        if log_file:
            print("\n\n..........................................................\n")
            print(log_message)
            log_file.write(log_message)
        return False
    
    return True


def remove_external_overhangs(overhangs, external_overhangs):
    """
    Remove external overhangs from a list of overhangs.

    We don't want to use these for assembly of inner fragments parts;
    so they are removed from overhang sets used as constraints used
    with goldenhinges.

    Parameters:
    overhangs (list): List of overhangs to filter.
    external_overhangs (list): List of external overhangs to remove.

    Returns:
    list: Filtered list of overhangs.
    """
    return [
        o for o in overhangs 
        if o not in external_overhangs and o not in (
            reverse_complement(oh) for oh in external_overhangs
        )
    ]


def create_chisels_df(
    accession, seq_id, original_sequence, chiseled_sequence, original_cai, chiseled_cai
):
    """
    Create a dataframe containing chiseled sequence information.

    Parameters:
    seq_id (str): Identifier for the sequence.
    original_sequence (str): Original DNA sequence.
    chiseled_sequence (str): Chiseled DNA sequence.
    original_cai (float): CAI for the original sequence.
    chiseled_cai (float): CAI for the chiseled sequence.

    Returns:
    pd.DataFrame: Dataframe containing the chiseled sequence information.
    """
    df = pd.DataFrame({
        'accession': [accession],
        'Seq_ID': [seq_id],
        'Original_sequence': [original_sequence],
        'Chiseled_sequence': [chiseled_sequence],
        'Original_CAI': original_cai,
        'Chiseled_CAI': chiseled_cai
    })
    return df

def check_membership(main_list, list_of_lists):
    """
    Check if the possible overhangs used for cloning are present in the 
    high fidelity overhang set being fed to goldenhinges as a constraint.
    """
    main_set = set(main_list)
    for sublist in list_of_lists:
        # Check if there's any overlap between the sublist and the main_set
        if not main_set.intersection(sublist):
            return False
    return True


def check_codon_usage_table(codon_usage_table):
    """
    Check the codon usage table for errors.

    Parameters:
    codon_usage_table (dict): Dictionary mapping amino acids to their codon 
                              frequencies.
    """
    for aa, codons in codon_usage_table.items():
        max_freq = max(codons.values())
        print(f"Max frequency for {aa}: {max_freq}")
        for codon, freq in codons.items():
            if freq > max_freq:
                print(
                    f"Error: Codon {codon} for {aa} has a frequency ({freq}) higher "
                    f"than the max frequency ({max_freq})"
                )


def calculate_cai(sequence, codon_usage_table):
    """
    Calculate the Codon Adaptation Index (CAI) for a given DNA sequence.

    The CAI is calculated by first breaking down the DNA sequence into its individual 
    codons. Each codon is then compared to a codon usage table, which provides the 
    frequency of each codon for its corresponding amino acid in the target organism.
    For each codon, a relative adaptiveness score is calculated by comparing the observed 
    codon frequency to the frequency of the most common codon for the same amino acid. 
    These relative adaptiveness scores are combined using a geometric mean, resulting in
    the CAI value, which ranges from 0 to 1. A higher CAI indicates that the sequence is 
    more optimally adapted to the organism's codon preferences.

    Parameters:
    sequence (str): DNA sequence.
    codon_usage_table (dict): Dictionary mapping amino acids to their codon 
                              frequencies.

    Returns:
    float: The CAI value, which is a measure of the sequence's adaptation to 
           the codon usage.

    See: Sharp, P. M., & Li, W. H. (1987). "The codon adaptation index—a measure 
    of directional synonymous codon usage bias, and its potential applications." 
    Nucleic Acids Research, 15(3), 1281-1295.
    
    """
    if len(sequence) % 3 != 0:
        raise ValueError("Sequence length should be a multiple of 3")

    # Split sequence into codons
    codons = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]

    # Prepare to collect log frequencies
    codon_log_frequencies = {}
    best_log_frequencies = {}

    # First pass: determine the best log frequencies
    for aa, codons_dict in codon_usage_table.items():
        max_freq = max(codons_dict.values(), default=1e-4)  # Avoid log(0) issue
        best_log_frequencies[aa] = np.log(max_freq)

    # Second pass: collect all log frequencies and calculate relative adaptiveness
    relative_adaptiveness = []
    for codon in codons:
        aa = next(
            (aa for aa, codons_dict in codon_usage_table.items() if codon in codons_dict), 
            None
        )
        if not aa:
            continue  # Skip if the codon doesn't correspond to any amino acid

        codon_freq = codon_usage_table[aa].get(codon, 1e-4)  # Again, prevent log(0)
        codon_log_frequencies[codon] = np.log(codon_freq)
        optimal_log_frequency = best_log_frequencies[aa]
        current_log_frequency = codon_log_frequencies[codon]
        ra = np.exp(current_log_frequency - optimal_log_frequency)
        relative_adaptiveness.append(ra)

    # Calculate the geometric mean of relative adaptiveness
    if relative_adaptiveness:
        geometric_mean = np.exp(np.mean(np.log(relative_adaptiveness)))
    else:
        geometric_mean = 0  # Handle case with no valid codons

    return geometric_mean

def calculate_fidelity_score(
    sequences: list, potapov_data: pd.DataFrame
) -> float:
    """
    Calculate the fidelity score for a list of overhangs based on Potapov 
    or other ligtation specificty data.

    The fidelity score calculation evaluates the specificity of a set of DNA 
    sequences by using on-target and off-target ligation data, such as that 
    from Potapov et al. For each sequence, the fidelity score is calculated by 
    comparing its off-target interactions with its on-target interactions. The 
    score for a sequence is determined as "1 minus the ratio of the off-target 
    score to the sum of the on-target and off-target scores." The final 
    fidelity score for a set of overhangs is the product of the individual 
    fidelity scores of each overhang. 

    see:
    Potapov, V., Ong, J. L., Kucera, R. B., Gushchina, I., Strychalski, E. A., 
    Schultz, D. N., & Potapov, V. (2018). Comprehensive profiling of four 
    base overhang ligation fidelity by T4 DNA Ligase and application to DNA 
    assembly. Nucleic Acids Research, 46(10), e79. doi:10.1093/nar/gky295

    Parameters:
    sequences (list): List of DNA sequences to be evaluated.
    potapov_data (pd.DataFrame): DataFrame containing Potapov data for on-target 
                                 and off-target scores.

    Returns:
    float: The final fidelity score for the provided sequences.

    """
    # Initialize a dictionary to store individual fidelity scores
    results = {"IndividualFidelity": []}

    # Filter out reverse complements and remove duplicates
    sequences = np.unique(filter_reverse_complements(sequences)).tolist()

    # Generate reverse complements for all unique sequences
    rc_sequences = [reverse_complement(seq) for seq in sequences]

    # Combine original and reverse complement sequences
    combined_sequences = list(set(sequences + rc_sequences))

    # Compute fidelity scores for each sequence/rc complement pair
    for seq in sequences:
        on_target_score, off_target_score = compute_target_scores(
            seq, combined_sequences, potapov_data
        )

        # Calculate individual fidelity as 1 - (off_target / (on_target + off_target))
        individual_fidelity = 1 - (
            off_target_score / (on_target_score + off_target_score)
        )

        # Append the individual fidelity score to results
        results["IndividualFidelity"].append(individual_fidelity)

    # Compute the final fidelity score as the product of all individual fidelities
    final_fidelity = np.prod(results["IndividualFidelity"])

    return final_fidelity


def compute_target_scores(
    seq: str, combined_sequences: list, potapov_data: pd.DataFrame
) -> tuple:
    """
    Compute the on-target and off-target scores for a given sequence based on 
    Potapov data.

    This function computes the on-target and off-target scores for a DNA 
    sequence based on Potapov data. It retrieves rows from the Potapov data 
    where the overhang matches either the sequence or its reverse complement. 
    The on-target score is determined by summing specific values from these 
    rows, while the off-target score is calculated by summing all values in 
    the relevant rows, subtracting the contributions from the on-target scores.
    The function returns these two scores as a tuple.

    Parameters:
    seq (str): The DNA sequence for which scores are to be computed.
    combined_sequences (list): List of combined sequences to be considered for 
                               scoring.
    potapov_data (pd.DataFrame): DataFrame containing Potapov data with overhang 
                                 sequences and scores.

    Returns:
    tuple: A tuple containing the on-target score and the off-target score.
    """

    # Compute the reverse complement of the given sequence.
    rc = reverse_complement(seq)

    # Remove duplicates from the combined sequences list.
    unique_combined_sequences = list(set(combined_sequences))

    # Extract rows from the Potapov data where the overhang matches the reverse 
    # complement of the sequence.
    rc_rows = potapov_data.loc[
        potapov_data["Overhang"] == rc, unique_combined_sequences
    ]

    # Extract rows from the Potapov data where the overhang matches the sequence.
    seq_rows = potapov_data.loc[
        potapov_data["Overhang"] == seq, unique_combined_sequences
    ]

    # Calculate the on-target score by summing the values from both the reverse 
    # complement and sequence rows.
    on_target_score = rc_rows[seq].values[0] + seq_rows[rc].values[0]

    # Calculate the off-target score by summing all values in the rows and 
    # subtracting the on-target contributions.
    off_target_score = (
        rc_rows.sum(axis=1).values[0] - rc_rows[seq].values[0] + 
        seq_rows.sum(axis=1).values[0] - seq_rows[rc].values[0]
    )

    return on_target_score, off_target_score


def reverse_complement(seq: str) -> str:
    """
    Compute the reverse complement of a DNA sequence.

    Parameters:
    seq (str): DNA sequence.

    Returns:
    str: Reverse complement of the DNA sequence.
    """
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement_dict[base] for base in reversed(seq))


def read_fasta(file_path):
    """
    Read sequences from a FASTA file.

    Parameters:
    file_path (str): Path to the FASTA file.

    Returns:
    list: List of sequences with their IDs.
    """
    file_path = f"{file_path}"
    sequences = []
    with open(file_path, 'r') as file:
        seq_id = None
        seq = ''
        for line in file:
            if line.startswith('>'):
                if seq_id:
                    sequences.append({'id': seq_id, 'sequence': seq})
                seq_id = line[1:].split(' ')[0].strip()
                seq = ''
            else:
                seq += line.strip()
        if seq_id:
            sequences.append({'id': seq_id, 'sequence': seq})
    return sequences


def filter_reverse_complements(sequences: List[str]) -> List[str]:
    """
    Filter out sequences that are reverse complements of each other, keeping 
    only one instance.

    Parameters:
    sequences (List[str]): List of DNA sequences to be filtered.

    Returns:
    List[str]: A list of sequences where reverse complements are filtered out.
    """
    # Remove duplicates while preserving order
    unique_sequences = []
    [unique_sequences.append(seq) for seq in sequences if seq not in unique_sequences]

    # Filter out sequences where a reverse complement is already present
    set1 = [
        seq for seq in unique_sequences if reverse_complement(seq) not in 
        unique_sequences or seq <= reverse_complement(seq)
    ]

    return set1

def print_organisms_from_file(filename):
    with open(filename, 'r') as f:
        for line in f:
            print(line.strip())

def write_sequences_to_fasta(sequences: pd.DataFrame, output_file: str):
    """
    Write sequences to a FASTA file.

    Parameters:
    sequences (pd.DataFrame): DataFrame containing the sequence IDs and 
                              sequences.
    output_file (str): Path to the output FASTA file.
    """
    with open(output_file, 'w') as file:
        for index, row in sequences.iterrows():
            seq_id = row['seq_id']
            number = row['number']
            sequence = row['chiseled_sequence']
            file.write(f">{seq_id}_chiseled_{number}\n")
            file.write(sequence + "\n\n")

def custom_warning_handler(message, category, filename, lineno, file=None, line=None):
    if category == BiopythonParserWarning:
        print(f"Warning: {message}")
        print("This warning indicates that there may be malformed lines in your GenBank file.")
        print("Please check the file and consider correcting any formatting issues.")
    else:
        warnings.showwarning(message, category, filename, lineno, file, line)


def annotate_genbank(input_file, output_file, yaml_file, id=None, pct=20):
    try:
        # Load the GenBank file
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            records = list(SeqIO.parse(input_file, "genbank"))
            for warning in w:
                custom_warning_handler(warning.message, warning.category, warning.filename, warning.lineno)

        if not records:
            raise ValueError("No records found in the GenBank file.")

        # Load the YAML file
        with open(yaml_file, 'r') as yf:
            yaml_data = yaml.safe_load(yf)

        # Extract species and codon_opt from YAML data
        species = yaml_data.get('species', 'none')
        codon_opt = yaml_data.get('codon_opt', 'none')

        # Annotate records using the parsed YAML data
        annotated_records = annotate_from_yaml(yaml_data, species, codon_opt, records, pct)

        # Write the modified records back to the output file
        with open(output_file, "w") as output_handle:
            SeqIO.write(annotated_records, output_handle, "genbank")

        # Post-process the output file to replace the specific label
        with open(output_file, "r") as file:
            content = file.read()

        # Perform the replacement
        content = content.replace(
            f'                     /label="~EnforceChanges(minimum_percent={pct})"',
            f'                     /label="@EnforceChanges(minimum={pct}%)"'
        )

        # Write the modified content back to the file
        with open(output_file, "w") as file:
            file.write(content)

        print(f"Annotation completed successfully. Output written to {output_file}")

    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        print("Please check the inputs and try again.")


def annotate_from_yaml(yaml_data, species, codon_opt, records, pct):
    try:
        regulatory_features = [
            "promoter", "terminator", "RBS", "regulatory",
            "rep_origin", "protein_bind", "tRNA", "oriT",
            "stem_loop", "enhancer", "5'UTR", "3'UTR",
            "polyA_signal", "LTR", "misc_binding", "tRNA",
            "rRNA", "misc_RNA"
        ]

        for record in records:
            new_features = []  # Collect new features here
            for feature in record.features:
                try:
                    if feature.type in regulatory_features:
                        avoid_changes_feature = SeqFeature(
                            feature.location,
                            type="misc_feature",
                            qualifiers={"label": "@AvoidChanges"}
                        )
                        new_features.append(avoid_changes_feature)

                    elif feature.type in ["CDS", "ORF"]:
                        misc_feature = SeqFeature(
                            feature.location,
                            type="misc_feature",
                            qualifiers={"label": "@EnforceTranslation"}
                        )
                        new_features.append(misc_feature)

                        if codon_opt == 'match_codon_usage':
                            codon_feature = SeqFeature(
                                feature.location,
                                type="misc_feature",
                                qualifiers={"label": f"~match_codon_usage(species={species})"}
                            )
                            new_features.append(codon_feature)

                        elif codon_opt == 'use_best_codon':
                            codon_feature = SeqFeature(
                                feature.location,
                                type="misc_feature",
                                qualifiers={"label": f"~use_best_codon(species={species})"}
                            )
                            new_features.append(codon_feature)

                        elif codon_opt == 'hybrid':
                            codon_feature = SeqFeature(
                                feature.location,
                                type="misc_feature",
                                qualifiers={"label": f"~use_best_codon(species={species})"}
                            )
                            new_features.append(codon_feature)
                            enforce_changes_feature = SeqFeature(
                                feature.location,
                                type="misc_feature",
                                qualifiers={"label": f"~EnforceChanges(minimum_percent={pct})"}
                            )
                            new_features.append(enforce_changes_feature)

                        elif codon_opt == 'harmonize_rca':
                            codon_feature = SeqFeature(
                                feature.location,
                                type="misc_feature",
                                qualifiers={"label": f"~match_codon_usage(species={species})"}
                            )
                            new_features.append(codon_feature)
                except Exception as e:
                    print(f"Error processing feature {feature.type}: {e}")
                    return records

            # Append all new features to the record
            record.features.extend(new_features)

            try:
                for constraint in yaml_data.get("constraints", []):
                    if constraint['type'] == "EnforceTranslation":
                        continue
                    label = create_feature_label("@", constraint)
                    misc_feature = SeqFeature(
                        FeatureLocation(0, len(record)),
                        type="misc_feature",
                        qualifiers={"label": label}
                    )
                    record.features.append(misc_feature)
            except Exception as e:
                print(f"Error processing constraints from YAML data: {e}")
                return records

            try:
                for optimization in yaml_data.get("optimizations", []):
                    label = create_feature_label("~", optimization)
                    misc_feature = SeqFeature(
                        FeatureLocation(0, len(record)),
                        type="misc_feature",
                        qualifiers={"label": label}
                    )
                    record.features.append(misc_feature)
            except Exception as e:
                print(f"Error processing optimizations from YAML data: {e}")
                return records

        return records
    except Exception as e:
        print(f"An error occurred during annotation: {e}")
        return records


def create_feature_label(prefix, data):
    try:
        label = f"{prefix}{data['type']}"
        additional_info = [
            f"{key}={value}" for key, value in data.items() if key != 'type'
        ]
        if additional_info:
            label += f"({', '.join(additional_info)})"
        return label
    except KeyError as e:
        print(f"Missing key in data: {e}")
        print("Please ensure the YAML data contains the correct keys.")
        return f"{prefix}unknown"

def adjust_feature_locations(features, offset):
    adjusted_features = []
    for feature in features:
        if feature.location is not None:
            start = feature.location.start + offset
            end = feature.location.end + offset
            new_location = FeatureLocation(
                start, end, strand=feature.location.strand
            )
            new_feature = SeqFeature(
                location=new_location, type=feature.type,
                qualifiers=feature.qualifiers
            )
            adjusted_features.append(new_feature)
    return adjusted_features


def log_and_print(message, log_file, quiet="off", width=80):
    lines = message.split('\n')
    wrapped_message = '\n'.join([
        textwrap.fill(line, width=width) for line in lines
    ])
    if quiet=="off":
        print(wrapped_message)

    log_file.write(wrapped_message + '\n')


def report_gb_cai(input_genbank, output_genbank, species):


    codon_usage_table = get_codons_table(species)

    def extract_cds_and_calculate_cai(genbank_file, codon_usage_table):
        records = SeqIO.parse(genbank_file, "genbank")
        data = []
        for record in records:
            for feature in record.features:
                if feature.type in ["CDS", "ORF"]:
                    seq = feature.location.extract(record).seq
                    cai = calculate_cai(str(seq), codon_usage_table)
                    annotation_fields = [
                        'standard_name', 'label', 'name', 'gene', 'id'
                    ]
                    annotation_value = None
                    for field in annotation_fields:
                        if field in feature.qualifiers:
                            annotation_value = feature.qualifiers[field][0]
                            break
                    data.append({
                        'ID': record.id,
                        'Annotation': annotation_value,
                        'Sequence': str(seq),
                        'CAI': cai
                    })
        return data

    input_data = extract_cds_and_calculate_cai(input_genbank, codon_usage_table)
    input_df = pd.DataFrame(input_data)

    output_data = extract_cds_and_calculate_cai(output_genbank, codon_usage_table)
    output_df = pd.DataFrame(output_data)

    return input_df, output_df


def save_cai_summary_to_file(input_df, output_df, output_path):
    with open(output_path, 'w') as f:
        f.write("Input GenBank CAI Summary:\n")
        f.write(input_df[['ID', 'CAI', 'Annotation']].to_string(index=False))
        f.write("\n\nOutput GenBank CAI Summary:\n")
        f.write(output_df[['ID', 'CAI', 'Annotation']].to_string(index=False))


def parse_and_load_config(default_yml):
    args = parse_arguments(default_yml=default_yml)
    default_values = set_defaults()
    updated_defaults, tag = load_config_and_set_globals(args, default_values)
    globals().update(vars(updated_defaults))
    return args, updated_defaults, tag


def log_and_exit(message, log_file, exit_code=1):
    log_and_print(message, log_file)
    log_file.write(f"Exit code: {exit_code}\n")
    log_file.write(
        f"Exit time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
    )
    sys.exit(exit_code)


def calculate_codon_frequencies(
        file_path='data/cleaned_coco.tsv', species_identifier=None):

    """

    DNAChisel's built-in codon table is based on the outdated Kazusa dataset and lacking
    data for many species and fails when the Kazusa website goes down. We use this as
    a function as a preferred alternative for generating codon usage tables. The function
    generates a usage table for a specified species using a trimmed down version of 
    the cocoputs dataset. The function loads codon frequency data by species name or NCBI 
    taxonomic ID. If a species name is provided, the function searches for it in the dataset 
    and selects the species with the highest number of codons if multiple matches are found; 
    this is becasue the cocoputs data set is large and often contains many entries for the
    same species. The codon frequencies are then calculated for each amino acid based on the 
    species-specific codon usage data. The result is a dictionary mapping amino acids to 
    their codon frequencies, along with the species name.

    See: Zhou, Z., Dang, Y., Zhou, M., Li, L., Yu, C. H., Fu, J., Chen, S., & Liu, Y. (2016). 
    Codon usage is an important determinant of gene expression levels largely through its 
    effects on transcription. Proceedings of the National Academy of Sciences, 113(41), 
    E6117-E6125.

    """

    codon_data = pd.read_csv(file_path, sep='\t')
    codon_data['Taxid'] = codon_data['Taxid'].astype(str)

    species_name = ""

    if isinstance(species_identifier, str) and not species_identifier.isdigit():
        matching_species = codon_data[
            codon_data['short_name'] == species_identifier]
        if len(matching_species) > 1:
            species_data = matching_species.loc[
                matching_species['# Codons'].idxmax()]
            print(
                f"Using Taxid={species_data['Taxid']} ({species_data['Species']}).")
            print(
                f"If {species_data['Species']} was not your intended species, "
                f"specify an NCBI taxid."
            )
            time.sleep(2)
        elif len(matching_species) == 0:
            print("No matching species found.")
            return None, None
        else:
            species_data = matching_species.iloc[0]
            print(
                f"Using Taxid={species_data['Taxid']} ({species_data['Species']}).")

    elif isinstance(species_identifier, (int, str)) and str(
            species_identifier).isdigit():
        species_data = codon_data[
            codon_data['Taxid'] == str(species_identifier)]
        if len(species_data) == 0:
            print("No matching species found.")
            return None, None
        else:
            species_data = species_data.iloc[0]
    else:
        print("Invalid species identifier.")
        return None, None

    species_name = f"{species_data['Species']} (taxid: {species_data['Taxid']})"
    codon_columns = codon_data.columns[8:]

    amino_acids = {
        '*': ['TAA', 'TAG', 'TGA'],
        'A': ['GCA', 'GCC', 'GCG', 'GCT'],
        'C': ['TGC', 'TGT'],
        'D': ['GAC', 'GAT'],
        'E': ['GAA', 'GAG'],
        'F': ['TTC', 'TTT'],
        'G': ['GGA', 'GGC', 'GGG', 'GGT'],
        'H': ['CAC', 'CAT'],
        'I': ['ATA', 'ATC', 'ATT'],
        'K': ['AAA', 'AAG'],
        'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
        'M': ['ATG'],
        'N': ['AAC', 'AAT'],
        'P': ['CCA', 'CCC', 'CCG', 'CCT'],
        'Q': ['CAA', 'CAG'],
        'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
        'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
        'T': ['ACA', 'ACC', 'ACG', 'ACT'],
        'V': ['GTA', 'GTC', 'GTG', 'GTT'],
        'W': ['TGG'],
        'Y': ['TAC', 'TAT']
    }

    codon_table = {}
    for aa, codons in amino_acids.items():
        aa_total = species_data[codons].sum()
        codon_table[aa] = {
            codon: round(species_data[codon] / aa_total, 3) for codon in codons
        }

    return codon_table, species_name


def write_fasta(df, filename):
    with open(filename, 'w') as f:
        for index, row in df.iterrows():
            f.write(f">{row['Seq_ID']}\n{row['Chiseled_sequence']}\n")


def read_fasta2(filename):
    try:
        records = list(SeqIO.parse(filename, "fasta"))
        sequences = {}
        for record in records:
            id_prefix = re.sub(r'.\d+$', '', record.id)
            if id_prefix not in sequences:
                sequences[id_prefix] = []
            sequences[id_prefix].append(str(record.seq))
        return sequences
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        return {}


def calculate_edit_distance(seq1, seq2):
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def calculate_distance_matrix(sequences):
    n = len(sequences)
    seq_length = len(sequences[0])
    distance_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            edit_distance = calculate_edit_distance(sequences[i], sequences[j])
            percent_distance = (edit_distance / seq_length) * 100
            distance_matrix[i, j] = percent_distance
            distance_matrix[j, i] = percent_distance

    return distance_matrix


def calculate_pairwise_identity(matrix):
    identity_matrix = 1 - (matrix / 100)
    return identity_matrix


def print_matrix(matrix, title, log_file, is_percentage=False, decimals=3):
    formatted_matrix = np.round(matrix, decimals)
    output = [title]
    if is_percentage:
        for row in formatted_matrix:
            output.append(
                ' '.join(f"{cell:.{decimals}f}%" for cell in row)
            )
    else:
        for row in formatted_matrix:
            output.append(
                ' '.join(f"{cell:.{decimals}f}" for cell in row)
            )

    for line in output:
        print(line)

    for line in output:
        log_file.write(line + '\n')


def calculate_average_distance(matrix):
    n = matrix.shape[0]
    np.fill_diagonal(matrix, np.nan)
    average_distance = np.nanmean(matrix)
    return average_distance


def call_tweakers_subsets(log_file, tweak_n, tweak_cai):
    print("\nTweaking selections initiated...")
    script_path = 'scripts/tweaker2.py'
    
    # Construct the command to call the script with the required arguments
    command = [
        'python', script_path,
        '--file', log_file,
        '--tweak_n', str(tweak_n),
        '--tweak_cai', str(tweak_cai)
    ]
    
    # Call the script
    result = subprocess.run(command, capture_output=True, text=True)
    


def call_tweaker_2(log_file, tweak_n):

    script_path = 'scripts/tweaker2.py'
    
    # Construct the command to call the script with the required arguments
    command = [
        'python', script_path,
        '--file', log_file,
        '--tweak_n', str(tweak_n),
    ]
    
    # Call the script
    result = subprocess.run(command, capture_output=True, text=True)
    
    # Print the output and error messages
    print(result.stdout)
    print(result.stderr)

def calculate_segment_length(pcr_5p_cut, primer_length, oligo_length):
    """
    Calculate the segment length based on the given parameters.

    Parameters:
    base_5p_end (str): The 5' base sequence.
    primer_length (int): The length of the primer.
    oligo_length (int): The total length of the oligo.

    Returns:
    int: The calculated segment length.
    """
    segment_length = oligo_length - (2 * len(pcr_5p_cut) + 2 * primer_length)
    return segment_length

