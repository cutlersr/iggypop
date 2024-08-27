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
import copy
import warnings
from pop_helpers import *
from headers import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import BiopythonParserWarning
from copy import deepcopy
from goldenhinges import OverhangsSelector
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
    EnforcePatternOccurence
)
from python_codon_tables import get_codons_table
from Bio.Align import PairwiseAligner

# Ignore warnings
warnings.simplefilter(action='ignore', category=Warning)

def chisel(
    sequence, method, codon_table, original_species, intron_constraints,
    loc_constraints_left, loc_constraints_right, left_bounds,
    right_bounds, log_file, reports, pct, last_try, repeat, quiet="off", 
    file="", results_path=""
):
    """
    Perform dnachisel optimizations using specified constraint and optimization
    parameters that are specified in a yaml config file or the command line;
    iggpop/yaml for config files.

    Parameters:
    sequence (str): DNA sequence to be optimized.
    method (str): Optimization method to be used (--codon_opt on the command line)
    codon_table (str): Codon usage table for the species.
    original_species (str): Original species for codon optimization (relevant for RCA).
    intron_constraints (list): List of intron junction sequence constraints.
    loc_constraints_left (list): List of left bounds for intron location constraints.
    loc_constraints_right (list): List of right bounds for intron location constraints.
    left_bounds (list): List of left bounds for fixing sequences at intron junctions.
    right_bounds (list): List of right bounds for fixing sequences at intron junctions.
    log_file (str): Path to the log file for logging messages.
    quiet (str): minimze prints to terminal.
    file (str): Path to the YAML file containing constraints and optimizations.
    results_path (str): Path to save the optimization report.

    Returns:
    str: Optimized DNA sequence.
    """
    file_path = f"{file}"  # Set the file path for the YAML file

    try:
        # Load data from the YAML file
        with open(file_path, 'r') as f:
            data = yaml.safe_load(f)
    except FileNotFoundError:
        log_and_print(f"YAML file not found: {file_path}", log_file, quiet)
        return None
    except yaml.YAMLError as e:
        log_and_print(f"Error parsing YAML file: {e}", log_file, quiet)
        return None

    try:
        print()
        # Parse constraints from the YAML file
        constraints = [
            eval(
                f"{entry['type']}("
                f"{', '.join(f'{key}={repr(value)}' for key, value in entry.items() if key != 'type')})"
            )
            for entry in data.get("constraints", [])
        ]
    except Exception as e:
        log_and_print(f"Error parsing constraints: {e}", log_file, quiet)
        return None

    try:
        # Parse optimizations from the YAML file
        optimizations = [
            eval(
                f"{entry['type']}("
                f"{', '.join(f'{key}={repr(value)}' for key, value in entry.items() if key != 'type')})"
            )
            for entry in data.get("optimizations", [])
        ]
    except Exception as e:
        log_and_print(f"Error parsing optimizations: {e}", log_file, quiet)
        return None

    # Hybrid mode uses use_best_codon w. EnforceChanges to create seq diversity
    # When running w/ deintronize mode, you only want it turned on for the first chisel
    if method == "hybrid" and last_try == 0:
        optimizations.append(
            CodonOptimize(
                codon_usage_table=codon_table, method='use_best_codon'
            )
        )
        constraints.append(
            EnforceChanges(minimum_percent=pct)
        )  # Enforce a minimum percent of changes

    elif method == "match_codon_usage" and repeat == 0:
        optimizations.append(
            CodonOptimize(
                codon_usage_table=codon_table,
                method=method,
                original_species=original_species
            )
        )

    elif method == "match_codon_usage" and repeat > 0:
        optimizations.append(
            CodonOptimize(
                codon_usage_table=codon_table,
                method=method,
                original_species=original_species
            )
        )
        constraints.append(
            EnforceChanges(minimum_percent=pct)
        )  # Enforce a minimum percent of changes

    elif method == "use_best_codon":
        optimizations.append(
            CodonOptimize(
                codon_usage_table=codon_table, method='use_best_codon'
            )
        )
    else:
        log_and_print("Skipping codon optimization", log_file, quiet)

    # Update optimizations to focus on potential intron areas
    if intron_constraints:
        optimizations = [
            AvoidPattern(
                pattern=intron_constraints[i],
                location=(loc_constraints_left[i], loc_constraints_right[i], 1)
            )
            for i in range(len(intron_constraints))
        ]

    try:
        # Create the DNA optimization problem
        problem = DnaOptimizationProblem(
            sequence=sequence,
            constraints=constraints,
            objectives=optimizations
        )
        log_and_print("\nBefore chisel optimizations:", log_file, quiet)
        log_and_print("----------------------------", log_file, quiet)
        log_and_print(problem.constraints_text_summary(), log_file, quiet)
        log_and_print(problem.objectives_text_summary(), log_file, quiet)

        problem.resolve_constraints(final_check=True)
        problem.max_random_iters = 1000
        try:

            if reports:
                try:
                    problem.optimize_with_report(target=results_path)
                except Exception as e:
                    log_and_print(
                        f"Issue with report generation; try running without reports: {e}",
                        log_file,
                        quiet
                    )
                    return None

            else:
                problem.optimize()

            log_and_print("After chisel optimizations:", log_file, quiet)
            log_and_print("---------------------------", log_file, quiet)
            log_and_print(problem.constraints_text_summary(), log_file, quiet)
            log_and_print(problem.objectives_text_summary(), log_file, quiet)
        except Exception as e:
            log_and_print(f"optimization failure: {e}", log_file, quiet)
            return None

    except Exception as e:
        log_and_print(f"Chiseling failed: {e}", log_file, quiet)
        return None

    return problem.sequence

def find_cut_solution(
    sequence, overhang_sets, radius, external_overhangs, segment_length,
    n_tries, potapov_data
):
    """
    Find the best cut solution for a given sequence based on fidelity scores.

    Parameters:
    sequence (str): The DNA sequence to be cut.
    overhang_sets (list): List of possible overhang sets.
    radius (int): Radius around the cuts.
    external_overhangs (list): List of external overhangs used for cloning.
    segment_length (int): Optimal length of each segment.
    n_tries (int): Total number of solutions to find.
    potapov_data (pd.DataFrame): Data for fidelity score calculations.

    Returns:
    tuple: Best solution, sequence sets, fidelity scores, counter value for the
           best solution.
    """

    x_tries = 50  # Maximum number of tries for each set size
    seq_sets = {}  # Dictionary to store sequence sets
    fidelity_scores = {}  # Dictionary to store fidelity scores
    overall_counter = 0  # Counter for overall tries
    set_size_counter = 0  # Counter for tries per set size
    last_set_size = None  # Variable to track the last set size

    # When hinging, break the input sequence into this many segments; the
    # calculation factors in the longest possible length (2 * radius). This
    # calculation ensures that no fragment will exceed the segment_length parameter
    segments = math.ceil(len(sequence) / (segment_length - radius * 2))

    # Extract all possible overhangs within the radii of the cuts
    seqs_around_cuts = extract_four_mers_around_cuts(sequence, segments, radius)

    for _, set_size, *overhangs in overhang_sets:
        # Reset set size counter if the set size changes
        if set_size != last_set_size:
            set_size_counter = 0
            last_set_size = set_size

        # Skip if the conditions for set size or overall tries are not met
        if set_size_counter >= x_tries or set_size < segments - 1:
            continue

        # Flatten the list of possible overhangs and remove external overhangs
        overhangs = [item for sublist in overhangs for item in sublist]
        overhangs = remove_external_overhangs(overhangs, external_overhangs)

        # Check to see if set of possible overhangs for the sequence being
        # hinged are present in the current overhangs read from the ohset file.
        # If not, increase the counter value; it will keep trying for x_tries
        # and then move onto the next set of overhangs
        if not check_membership(overhangs, seqs_around_cuts):
            set_size_counter += 1
            continue

        try:
            selector = OverhangsSelector(
                possible_overhangs=overhangs, time_limit=100
            )
            solutions = selector.cut_sequence(
                sequence, equal_segments=segments, max_radius=radius,
                include_extremities=False, allow_edits=False, solutions=1,
                optimize_score=False
            )

            if solutions:
                sequences = [sol['sequence'] for sol in solutions]
                sequences.extend(external_overhangs)
                seq_sets[overall_counter] = {
                    'solutions': solutions, 'counter': overall_counter
                }
                fidelity_scores[overall_counter] = calculate_fidelity_score(
                    sequences, potapov_data
                )
                overall_counter += 1

        except Exception as e:
            print(f"Error during sequence cutting: {e}")
            continue

        set_size_counter += 1

        # Break the loop if the overall counter exceeds or equals n_tries
        if overall_counter >= n_tries:
            break

    if not seq_sets:
        print("No solutions found")
        return None, None, None, overall_counter

    best_fidelity_score = max(fidelity_scores.values(), default=0)
    best_solution_index = max(fidelity_scores, key=fidelity_scores.get, default=-1)
    best_solution = seq_sets.get(best_solution_index, {}).get('solutions')
    counter_value_for_best_solution = seq_sets.get(
        best_solution_index, {}).get('counter', 0
    )

    return best_solution, seq_sets, fidelity_scores, counter_value_for_best_solution


def get_overhang_sets(filename, external_overhangs):
    """
    Load overhang sets from the high fidelity oh_sets file and filter them to 
    remove the external overhangs used for cloning.

    Parameters:
    filename (str): Path to the Excel file containing overhang sets.
    external_overhangs (list): List of external overhangs to be filtered out.

    Returns:
    list: List of tuples containing external overhangs, set size, and 
          filtered overhang set.
    """
    # Load the Excel file into a pandas DataFrame
    df = pd.read_excel(filename)

    # Initialize an empty list to store overhang sets
    overhang_sets = []

    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        # Extract values from the current row
        external_overhangs = row['external_overhangs']
        set_size = row['set_size']
        overhang_set = row['hf_oh_set'].split(', ')

        # Filter out overhangs that are in the external_overhangs list
        filtered_overhang_set = [
            oh for oh in overhang_set if oh not in external_overhangs
        ]

        # Add the tuple (external_overhangs, set_size, filtered_overhang_set) 
        # to the list
        overhang_sets.append((external_overhangs, set_size, 
                              filtered_overhang_set))

    return overhang_sets

# Other utility functions

def extract_four_mers_around_cuts(sequence, segments, radius):
    """
    Extract four-mers around cut points in the sequence.

    Parameters:
    sequence (str): The DNA sequence.
    segments (int): Number of segments to divide the sequence into.
    radius (int): Radius around the cut points.

    Returns:
    list: List of four-mers around each cut point.
    """
    bp = 4  # base pairs in overhang
    segment_length = len(sequence) // segments
    cut_points = [
        (i * segment_length, (i + 1) * segment_length) for i in range(segments - 1)
    ]
    sequences_around_cuts = []
    for start, end in cut_points:
        segment = sequence[max(0, end - radius):min(len(sequence), end + radius)]
        sequences_around_cuts.append(segment)
    result = []
    for segment in sequences_around_cuts:
        xmers = [segment[i:i + bp] for i in range(len(segment) - bp + 1)]
        result.append(xmers)
    return result

def remove_external_overhangs(overhangs, external_overhangs):
    """
    Remove external overhangs from a list of overhangs.

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

def calculate_fidelity_score(
    sequences: list, potapov_data: pd.DataFrame
) -> float:
    """
    Calculate the fidelity score for a list of overhangs based on Potapov 
    data.

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

def filter_reverse_complements(sequences: List[str]) -> List[str]:
    """
    Filter out sequences that are reverse complements of each other, keeping 
    only one instance 
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

def compute_target_scores(
    seq: str, combined_sequences: list, potapov_data: pd.DataFrame
) -> tuple:
    """
    Compute the on-target and off-target scores for a given sequence based on 
    Potapov data.

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

