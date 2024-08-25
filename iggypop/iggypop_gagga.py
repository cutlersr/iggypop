import argparse
import itertools
import numpy as np
import pandas as pd
import random
import datetime
import os
from deap import base, creator, tools, algorithms
from typing import Tuple, List
from collections import OrderedDict


# parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description='GAGGA: GA-based Golden Gate Assembly overhang set optimizer')
    parser.add_argument('--potapov_data', type=str, default="data/FileS03_T4_18h_25C.xlsx", 
                        help='Potapov data file to use for calculating fidelity, on-target, and off-target hybridization scores; other options are present in the data folder')
    parser.add_argument('--set_size', type=int, default=30, help='Overhang set size; default = 30')
    parser.add_argument('--n_best', type=int, default=50, help='Number of sets to report in output; default = 50')
    parser.add_argument('--fixed_overhangs', type=str, nargs='+', default=["AATG", "GCTT"],
                        help='Fixed overhangs; these will be present in every set generated')
                        # These are the lowest scoring hybridization sequences and the most promiscuous cross-binding overhangs
    parser.add_argument('--forbidden_sequences', type=str, nargs='+', 
                        default=['GGAC', 'CGCG', 'GGCC', 'GCCA', 'CCGC', 'GGCG', 'GTCG', 'TCGC', 'GGGT', 'GGGG', 'TGCG', 
                        'TAAA', 'GCAT', 'GGCT', 'TTTA', 'GCGG', 'TATA', 'GGAT', 'CCCC', 'GGGC', 'GCGT', 'TTAA', 'GGAG', 'GGTG'],                    
                        help='Sequences to omit from solution sets; default = the 10 worst seqs in Potapov dataset')

    parser.add_argument('--pop_size', type=int, default=750, help='Population size for GA; default 500; reduce to speed up')
    parser.add_argument('--ngen', type=int, default=100, help='Number of generations')
    parser.add_argument('--mutpb', type=float, default=0.2, help='Mutation probability')
    parser.add_argument('--cxpb', type=float, default=0.4, help='Crossover probability')
    parser.add_argument('--alpha', type=float, default=2, help='Penalty parameter; higher values are more forigiving of small set sizes; increase for larger sets; default 2')
    parser.add_argument('--beta', type=float, default=2, help='Penalty parameter; higher values are more forigiving of small set sizes; increase for larger sets default 2')
    parser.add_argument('--indpb', type=float, default=0.05, help='Prob. of an overhang w/in a set being mutated; default = 0.05')
    parser.add_argument('--treshold', type=float, default=0.6, help='Similarlity threshold for increasing mutation/crossover rates; default = 0.6')
    parser.add_argument('--topn', type=int, default=50, help='Number of top individuals to include in set-similarity calculation; default = 50')
    parser.add_argument('--max_cxpb', type=float, default=0.8, help='Maxmimum crossover rate; default = 0.80')
    parser.add_argument('--max_mutpb', type=float, default=0.75, help='Maxmimum mutation rate; default = 0.875')
    parser.add_argument('--max_indpb', type=float, default=.075, help='Maxmimum indpb value ; default = 0.1')
    parser.add_argument('--min_improve', type=float, default=.0025, help='Stop run if minimum average improvement over last tries is < this value ; default = .01')
    parser.add_argument("--stag", type=int, default=10, help="Stop if score does not exceed 'min_improve' after this many generations; default = 10")
    parser.add_argument('--elite_set', type=str, default='none', 
                        help='Use one of the Potapov predefined sets to start the GAGGA-ing; options: 15, 20, 25, 30, user, none (default)')
    parser.add_argument('--user_elite_set', type=str, default='20', 
                        help='Comma-separated list of starting sequences to use as initiating seqs')
    parser.add_argument("--tournament_size", type=int, default=2, help="Tournament size for selection; default = 10")
 

    return parser.parse_args()

# Parse the arguments
args = parse_arguments()

# CONSTANTS & PARAMETERS
POTAPOV_DATA = args.potapov_data
SET_SIZE = args.set_size
N_BEST = args.n_best
FIXED_OVERHANGS = args.fixed_overhangs
FORBIDDEN_SEQUENCES = args.forbidden_sequences
POP_SIZE = args.pop_size
NGEN = args.ngen
MUTPB = args.mutpb
CXPB = args.cxpb
ALPHA = args.alpha
BETA = args.beta
TRESHOLD = args.treshold
INDPB = args.indpb
MAX_CXPB = args.max_cxpb
MAX_MUTPB = args.max_mutpb
MAX_INDPB = args.max_indpb
TOPN = args.topn
MIN_IMPROVEMENT = args.min_improve
ELITE_SET = args.elite_set
USER_ELITE_SET = args.user_elite_set
n_stagnant_generations = args.stag


# HELPER FUNCTIONS
def reverse_complement(seq: str) -> str:
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement_dict[base] for base in reversed(seq))


def compute_target_scores(seq: str, combined_sequences: list, potapov_data: pd.DataFrame) -> tuple:
    """
    Calculate on-target and off-target scores based on the given sequence and data.
    
    Parameters:
    - seq (str): The sequence for which the scores are to be calculated.
    - combined_sequences (list): A list of sequences.
    - potapov_data (pd.DataFrame): A dataframe containing overhang data and scores.
    
    Returns:
    - tuple: A tuple containing on-target and off-target scores.
    """
    
    # Get the reverse complement of the sequence
    rc = reverse_complement(seq)
    
    # Ensure the combined list of sequences and reverse complements is unique
    unique_combined_sequences = list(set(combined_sequences + [reverse_complement(seq) for seq in combined_sequences]))
    
    # Store the rows we're interested in to avoid repeated lookups
    rc_rows = potapov_data[potapov_data["Overhang"] == rc]
    seq_rows = potapov_data[potapov_data["Overhang"] == seq]
    
    # Calculate on_target_score directly using the stored rows
    on_target_score = rc_rows[seq].values[0] + seq_rows[rc].values[0]
    
    # Calculate off_target_score only for the pairings with the reverse complement of the sequence
    off_target_score = rc_rows[unique_combined_sequences].drop(columns=[seq]).sum().sum() + seq_rows[unique_combined_sequences].drop(columns=[rc]).sum().sum()
    
    return on_target_score, off_target_score

def reverse_complement(seq: str) -> str:
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement_dict[base] for base in reversed(seq))

def calculate_redundancy(sequences: list) -> tuple:
    # we want to penalize sets that have repeated sequences or palidromes
    # this function retyurns the number of unique non-palidromic / non-self pairing sequences in the set
    # sequences            
    unique_set = set()
    #added this to try an fix the omission issue
    sequences = list(sequences)
    for seq in sequences:
        rc = reverse_complement(seq)
        if seq == rc:  # Sequence is palindromic
            unique_set.add(seq)
        else:
            if seq not in unique_set and rc not in unique_set:
                unique_set.add(seq)
    #a set has redundant members when there are repeated sequences, reverse complement pairs, or palidromes            
    redundant_count = len(sequences) - len(unique_set) + sum(1 for seq in unique_set if seq == reverse_complement(seq))
    normalized_redundancy = redundant_count / len(sequences)
    return redundant_count, normalized_redundancy
  
  
def calculate_fidelity_score(sequences: list, potapov_data: pd.DataFrame) -> float:
    results = {"Sequence": [], "OnTarget": [], "OffTarget": [], "IndividualFidelity": []}
    #make sure that the fixed_overhangs are scored (not sure why but sometimes they are being left out)
    sequences = list(sequences) + FIXED_OVERHANGS
    sequences = np.unique(filter_reverse_complements(sequences)[0]).tolist()
    rc_sequences = [reverse_complement(seq) for seq in sequences]
    combined_sequences = list(set(sequences + rc_sequences))
    for seq in sequences:
        on_target_score, off_target_score = compute_target_scores(seq, combined_sequences, potapov_data)
        individual_fidelity = 1 - (off_target_score / (on_target_score + off_target_score))
        results["Sequence"].append(seq)
        results["OnTarget"].append(on_target_score)
        results["OffTarget"].append(off_target_score)
        results["IndividualFidelity"].append(individual_fidelity)
    final_fidelity = np.prod(results["IndividualFidelity"])
    return final_fidelity


def calculate_penalized_fidelity_score(sequences: list, potapov_data: pd.DataFrame, alpha: float, beta: float) -> float:
    raw_fidelity = calculate_fidelity_score(sequences, potapov_data)
    _, normalized_redundancy = calculate_redundancy(sequences)
    penalty_factor = np.exp(-alpha * (normalized_redundancy ** beta))
    penalized_fidelity = raw_fidelity * penalty_factor
    return penalized_fidelity


def filter_reverse_complements(sequences: List[str]) -> Tuple[List[str], List[str], List[str], List[str]]:
    def reverse_complement(seq: str) -> str:
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return ''.join([complement[base] for base in reversed(seq)])

    # Remove duplicates while preserving order
    sequences = list(sequences)

    unique_sequences = []
    [unique_sequences.append(seq) for seq in sequences if seq not in unique_sequences]

    set1 = [seq for seq in unique_sequences if reverse_complement(seq) not in unique_sequences or seq <= reverse_complement(seq)]
    set2 = [seq for seq in unique_sequences if reverse_complement(seq) not in unique_sequences or seq > reverse_complement(seq)]

    set1_no_palindromes = [seq for seq in set1 if seq != reverse_complement(seq)]
    set2_no_palindromes = [seq for seq in set2 if seq != reverse_complement(seq)]

    return set1, set2, set1_no_palindromes, set2_no_palindromes



# GENETIC ALGORITHM FUNCTIONS
def init_individual(set_size: int) -> list:
    individual = FIXED_OVERHANGS.copy() 
    # Remove FIXED_OVERHANGS and their reverse complements from possible_overhangs 
    rc_fixed = [reverse_complement(seq) for seq in FIXED_OVERHANGS]
    combined_fixed_rc = FIXED_OVERHANGS + rc_fixed
    available_overhangs = [oh for oh in possible_overhangs if oh not in combined_fixed_rc]
    individual += random.sample(available_overhangs, set_size - len(FIXED_OVERHANGS))
    return creator.Individual(individual)

def custom_mutation(individual: list, indpb: float) -> tuple: 
    # Remove FIXED_OVERHANGS from possible_overhangs to avoid duplicates
    rc_fixed = [reverse_complement(seq) for seq in FIXED_OVERHANGS]
    combined_fixed_rc = FIXED_OVERHANGS + rc_fixed
    available_overhangs = [oh for oh in possible_overhangs if oh not in combined_fixed_rc]
    #sample new seqs using the indpb parameter 
    for i in range(len(FIXED_OVERHANGS), len(individual)):
        if random.random() < indpb:
            individual[i] = random.choice(available_overhangs)
    return individual,


def custom_crossover(ind1: list, ind2: list) -> tuple:
    # Remove FIXED_OVERHANGS from possible_overhangs to avoid duplicates
    rc_fixed = [reverse_complement(seq) for seq in FIXED_OVERHANGS]
    combined_fixed_rc = FIXED_OVERHANGS + rc_fixed
    available_overhangs = [oh for oh in possible_overhangs if oh not in combined_fixed_rc]
    size = min(len(ind1), len(ind2))
    cxpoint1, cxpoint2 = sorted(random.sample(range(len(FIXED_OVERHANGS), size), 2))
    ind1[cxpoint1:cxpoint2], ind2[cxpoint1:cxpoint2] = ind2[cxpoint1:cxpoint2].copy(), ind1[cxpoint1:cxpoint2].copy()
    for idx, val in enumerate(ind1):
        if val in forbidden_sequences:
            ind1[idx] = random.choice(available_overhangs)
    for idx, val in enumerate(ind2):
        if val in forbidden_sequences:
            ind2[idx] = random.choice(available_overhangs)
    return ind1, ind2


def evaluate(individual: list) -> tuple:
    sequences = individual
    return calculate_penalized_fidelity_score(sequences, df, ALPHA, BETA),


def adjust_elite_set(elite_set: List[str], fixed_overhangs: List[str], forbidden_sequences: List[str], target_size: int) -> List[str]:
    # Remove fixed_overhangs and forbidden_sequences from elite_set
    rc_fixed = [reverse_complement(seq) for seq in FIXED_OVERHANGS]
    combined_fixed_rc = FIXED_OVERHANGS + rc_fixed
    adjusted_elite_set = [seq for seq in elite_set if seq not in combined_fixed_rc + FORBIDDEN_SEQUENCES]
    
    # Generate a list of allowable sequences
    all_sequences = ["".join(seq) for seq in itertools.product("ACGT", repeat=4)]
    allowable_sequences = [seq for seq in all_sequences if seq not in combined_fixed_rc + FORBIDDEN_SEQUENCES + adjusted_elite_set]
    
    # Adjust the size of the adjusted_elite_set
    while len(adjusted_elite_set) + len(fixed_overhangs) < target_size:
        # If adjusted_elite_set is smaller than target_size, add random sequences from allowable_sequences
        adjusted_elite_set.append(random.choice(allowable_sequences))
        allowable_sequences.remove(adjusted_elite_set[-1])
    while len(adjusted_elite_set) + len(fixed_overhangs) > target_size:
        # If adjusted_elite_set is larger than target_size, randomly remove sequences
        removed_seq = random.choice(adjusted_elite_set)
        adjusted_elite_set.remove(removed_seq)
    
    # Add fixed_overhangs to the beginning of adjusted_elite_set
    adjusted_elite_set = fixed_overhangs + adjusted_elite_set
    
    return adjusted_elite_set


def process_and_save_top_solutions_to_excel(top_solutions: list, potapov_data: pd.DataFrame, output_file: str) -> None:
    # Convert top_solutions to a set of unique individuals and then convert it back to a list
    top_solutions = list({tuple(ind) for ind in top_solutions})
    
    # Extract fidelity scores and non-redundant sequence counts for the unique individuals
    penalized_scores = [calculate_penalized_fidelity_score(ind, potapov_data, ALPHA, BETA) for ind in top_solutions]
    raw_scores = [calculate_fidelity_score(ind, potapov_data) for ind in top_solutions]
    penalties = [1 - (penalized / raw) if raw != 0 else 0 for penalized, raw in zip(penalized_scores, raw_scores)]
    
    # Sort unique individuals by raw scores in descending order
    sorted_indices = np.argsort(raw_scores)[::-1]
    sorted_solutions = [top_solutions[i] for i in sorted_indices]
    
    # Select the top 50 unique individuals with the highest raw scores
    sorted_solutions = sorted_solutions[:50]

    # Recalculate penalized_scores, raw_scores, and penalties for the top 25 individuals
    penalized_scores = [calculate_penalized_fidelity_score(ind, potapov_data, ALPHA, BETA) for ind in sorted_solutions]
    raw_scores = [calculate_fidelity_score(ind, potapov_data) for ind in sorted_solutions]
    penalties = [1 - (penalized / raw) if raw != 0 else 0 for penalized, raw in zip(penalized_scores, raw_scores)]
    
    # Create DataFrames for each set of sequences
    set1_df, set2_df, set1_no_palindromes_df, set2_no_palindromes_df = [], [], [], []
    for ind in sorted_solutions:
        unique_ind = list(set(ind[2:]))
        ext_oh = list(set(ind[0:2]))
        set1, set2, set1_no_palindromes, set2_no_palindromes = filter_reverse_complements(unique_ind)
        set1 = ext_oh + set1
        set1_df.append(', '.join(np.unique(set1)))
        set2_df.append(', '.join(np.unique(set2)))
        set1_no_palindromes_df.append(', '.join(np.unique(set1_no_palindromes)))
        set2_no_palindromes_df.append(', '.join(np.unique(set2_no_palindromes)))

    # Calculate similarity to top-scoring set
    top_scoring_set = set(sorted_solutions[0])
    similarities = [len(top_scoring_set.intersection(set(ind))) / SET_SIZE for ind in sorted_solutions]

    # Create a DataFrame to store the results

    from collections import OrderedDict
    results_df = pd.DataFrame({
        'hf_oh_set': set1_df,
        'fidelity': ['{:.1%}'.format(score) for score in raw_scores]
    })
    # Remove duplicates based on the 'hf_oh_set' column, keeping the first occurrence.
    results_df = results_df.drop_duplicates(subset='hf_oh_set', keep='first')


    # Create a DataFrame for the parameters
    parameters_data = {
        'SET_SIZE': SET_SIZE,
        'POP_SIZE': POP_SIZE,
        'TOURNSIZE': args.tournament_size,
        'MIN_IMPROVEMENT': MIN_IMPROVEMENT,
        'ALPHA': ALPHA,
        'BETA': BETA,
        'NGEN': NGEN,
        'MUTPB': MUTPB,
        'CXPB': CXPB,
        'INDPB': INDPB, 
        'MAX_CXPB': MAX_CXPB,
        'MAX_MUTPB': MAX_MUTPB,
        'MAX_INDPB': MAX_INDPB,
        'ELITE_SET': ELITE_SET,
        'USER_ELITE_SET': USER_ELITE_SET,
        'FIXED_OVERHANGS': ', '.join(FIXED_OVERHANGS),
        'FORBIDDEN_SEQUENCES': ', '.join(FORBIDDEN_SEQUENCES),
        'N_BEST': N_BEST,
        'TRESHOLD': TRESHOLD,
        'TOPN': TOPN, 
        'POTAPOV_DATA': POTAPOV_DATA,
        'GAGGA_SCRIPT': os.path.basename(__file__),
    }
    parameters_df = pd.DataFrame.from_dict(parameters_data, orient='index', columns=['Value']).T

    # Generate a timestamp
    timestamp = datetime.datetime.now().strftime('%H%M%S')

    # Modify the output filename to include set_size and timestamp
    filename_base, filename_extension = os.path.splitext(output_file)
    version = 1

    while os.path.exists(output_file):
        output_file = f"out/gagga/{filename_base}_{SET_SIZE}_{FIXED_OVERHANGS}_{NGEN}_gens_{timestamp}_v{version}{filename_extension}"
        version += 1    
    print(output_file)

    # Save the DataFrames to an Excel file
    with pd.ExcelWriter(output_file) as writer:
        results_df.to_excel(writer, sheet_name='Results', index=False)
        parameters_df.to_excel(writer, sheet_name='Parameters', index=False)

def is_non_redundant(individual, hof_items):
    """
    Check if the individual is non-redundant in the hall of fame.
    
    Parameters:
    - individual: The individual to check.
    - hof_items: A list of individuals already in the hall of fame.
    
    Returns:
    - True if the individual is non-redundant, False otherwise.
    """
    ind_set = set(individual)  # Convert the individual to a set
    for hof_ind in hof_items:
        hof_ind_set = set(hof_ind)  # Convert each hall of fame individual to a set
        if ind_set == hof_ind_set:
            return False  # Found a redundant individual
    return True  # No redundancy found

def eaSimpleWithElitism(population, toolbox, cxpb, mutpb, indpb, ngen, stats=None, halloffame=None, verbose=__debug__):
    logbook = tools.Logbook()
    logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

    # Define early stopping parameters
    max_scores = []

    # Evaluate the fitness of each individual in the population
    fitnesses = toolbox.map(toolbox.evaluate, population)
    for ind, fit in zip(population, fitnesses):
        ind.fitness.values = fit

    if halloffame is None:
        raise ValueError("halloffame parameter must not be empty!")

    halloffame.update(population)
    hof_size = len(halloffame.items) if halloffame else 0

    record = stats.compile(population) if stats else {}
    logbook.record(gen=0, nevals=len(population), **record)
    if verbose:
        print(logbook.stream)

    # Begin the generational process
    for gen in range(1, ngen + 1):
        # Select the next generation individuals
        offspring = toolbox.select(population, len(population) - hof_size)

        # Vary the pool of individuals with crossover and mutation
        offspring = algorithms.varAnd(offspring, toolbox, cxpb, mutpb)

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # Update the hall of fame with a check for redundancy
        new_hof_items = []
        for ind in offspring:
            if is_non_redundant(ind, halloffame.items):  # You need to implement this method
                new_hof_items.append(ind)
        halloffame.update(new_hof_items)

        # Check for early stopping
        max_scores.append(max(ind.fitness.values[0] for ind in population))
        if len(max_scores) > n_stagnant_generations:
            max_scores.pop(0)
            improvement = (max_scores[-1] - max_scores[0]) / max_scores[0]
            if improvement < MIN_IMPROVEMENT:
                break

        # Replace the current population by the offspring
        population[:len(offspring)] = offspring
        population[len(offspring):] = halloffame.items

        # Append the current generation statistics to the logbook
        record = stats.compile(population) if stats else {}
        logbook.record(gen=gen, nevals=len(invalid_ind), **record)
        if verbose:
            print(logbook.stream)

        # Get the highest-scoring individual
        highest_scoring_individual = halloffame.items[0]


        # Calculate the average similarity of the top N individuals
        topn_individuals = tools.selBest(population, k=TOPN)
        total_similarity = 0
        for i in range(len(topn_individuals)):
            for j in range(i+1, len(topn_individuals)):
                common_sequences = set(topn_individuals[i]) & set(topn_individuals[j])
                total_similarity += len(common_sequences) / SET_SIZE
        average_similarity = total_similarity / (len(topn_individuals) * (len(topn_individuals) - 1) / 2)
#        print(f"Average_similarity of the top {TOPN} individuals:, {average_similarity:.3f}")

        
        # Increase mutation and crossover rates if average similarity is above a threshold
        # if average_similarity > TRESHOLD:
        #     cxpb = min(MAX_CXPB, cxpb * 1.04)
        #     mutpb = min(MAX_MUTPB, mutpb * 1.03)
        #     indpb = min(MAX_INDPB, indpb * 1.05)


        # # Increase mutation and crossover rates if average similarity is above a threshold
        # if average_similarity > TRESHOLD:
        #     cxpb = min(MAX_CXPB, cxpb * 1.03)
        #     mutpb = min(MAX_MUTPB, mutpb * 1.02)
        #     indpb = min(MAX_INDPB, indpb * 1.04)

        # Increase mutation and crossover rates if average similarity is above a threshold
        if average_similarity > TRESHOLD:
            cxpb = min(MAX_CXPB, cxpb * 1.02)
            mutpb = min(MAX_MUTPB, mutpb * 1.02)
            indpb = min(MAX_INDPB, indpb * 1.03)



        # Print mutpb and cxpb values
#        print("mutpb:", f'{mutpb:.2f}', "cxpb:", f'{cxpb:.2f}')
#        print("")

    return


# MAIN 
df = pd.read_excel(POTAPOV_DATA)

forbidden_sequences = FORBIDDEN_SEQUENCES
possible_overhangs = [seq for seq in df.columns[1:] if seq not in forbidden_sequences]

creator.create("FitnessMulti", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMulti)

toolbox = base.Toolbox()

toolbox.register("individual", init_individual, SET_SIZE)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("mate", custom_crossover)
toolbox.register("mutate", custom_mutation, indpb=INDPB)
toolbox.register("select", tools.selTournament, tournsize=args.tournament_size)
toolbox.register("evaluate", evaluate)

# Parse tournament_size argument; if not provided, set it to pop_size/10
if args.tournament_size is None:
    args.tournament_size = args.pop_size // 10

# Parse elite_set argument
elite_set = ELITE_SET
if elite_set == 'user':
    # If elite_set is 'user', get the elite set from the user_elite_set argument
    elite_set = USER_ELITE_SET.split(',')
elif elite_set in ['15', '20', '25', '30']:
    # If elite_set is '15', '20', '25', or '30', get the elite set from the predefined sets
    predefined_elite_sets = {
        '15': ["TGCC", "GCAA", "ACTA", "TTAC", "CAGA", "TGTG", "GAGC", "CGAA",
           "AGGA", "ATTC", "ATAG", "AAGG", "AACT", "AAAA", "ACCG"],  

        '20': ["AGTG", "CAGG", "ACTC", "AAAA", "AGAC", "ATAG",
           "AACC", "TACA", "TAGA", "ATGC", "GATA", "GTAA", "CTGA",
           "ACAA", "AGGA", "ATTA", "ACCG", "GCGA", "CGAA", "CTCC"], 

        '25': ["CCTC", "CTAA", "GACA", "GCAC", "AATC", "GTAA", "TGAA",
           "ATTA", "AGGA", "ACAA", "TAGA", "CGGA", "CATA", "CAGC",
           "AACG", "AAGT", "AGAT", "ACCA", "AGTG", "GGTA", "GCGA",
           "AAAA", "ATGA", "CTCC", "CCAG"], 

        '30': ["CCAG", "AGAG", "TACA", "CTAA", "GGAA", "GCCA", "ACTC", "CTTC",
           "TCAA", "GATA", "ACTG", "AACT", "AAGC", "CATA", "GACC", "AGGA",
           "ATCG", "ATTA", "CGGA", "TAGA", "AGCA", "TGAA", "ACAT",
           "GTGA", "ACGA", "ATAC", "AAAA", "AAGG", "CAAC"],  
    }
    elite_set = predefined_elite_sets[elite_set]
else:
    elite_set = []


# Adjust the elite set
adjusted_elite_set = adjust_elite_set(elite_set, FIXED_OVERHANGS, FORBIDDEN_SEQUENCES, SET_SIZE)

# Create the initial population
population = toolbox.population(n=POP_SIZE)

# If adjusted_elite_set is not empty, set the first individual and generate the rest
if adjusted_elite_set:
    toolbox.register("individual", init_individual, len(adjusted_elite_set))
    
    population[0][:] = adjusted_elite_set
    for ind in population[1:]:
        ind[:] = toolbox.individual()
else:
    # If adjusted_elite_set is empty, register the "individual" function and generate all individuals randomly
    toolbox.register("individual", init_individual, SET_SIZE)
    for ind in population:
        ind[:] = toolbox.individual()

stats = tools.Statistics(lambda ind: ind.fitness.values)
stats.register("avg", np.mean, axis=0)
stats.register("std", np.std, axis=0)
stats.register("min", np.min, axis=0)
stats.register("max", np.max, axis=0)

halloffame = tools.HallOfFame(maxsize=50)


eaSimpleWithElitism(population, toolbox, CXPB, MUTPB, INDPB, NGEN, stats=stats, halloffame=halloffame)

# Retrieving the top solutions after the genetic algorithm run
top_solutions = tools.selBest(population, k=N_BEST)

# Create the results directory if it doesn't exist
if not os.path.exists(f'out/gagga'):
    os.makedirs(f'out/gagga')

# Post-processing and saving results
process_and_save_top_solutions_to_excel(top_solutions, df, 'gagga.xlsx')


