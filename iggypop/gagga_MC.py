#!/usr/bin/env python3

# Monte Carlo–based optimizer to select high‑affinity DNA overhang sets
# under user‑specified fixed and forbidden constraints.
# The default set of forbidden overhangs include all 4-mer pallindromes
# and the 10overhangs with the lowest fidelity scores in Potapov et al.'s data set

import pandas as pd
import numpy as np
import random
import math
import os
from openpyxl import load_workbook
import argparse
from datetime import datetime

# ----------------------------------------
# Utility functions
# ----------------------------------------

def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complement[base] for base in reversed(sequence))

def generate_all_4mers():
    """Generate all possible 4‑base DNA sequences (4‑mers)."""
    bases = ['A', 'T', 'C', 'G']
    return [a+b+c+d for a in bases for b in bases for c in bases for d in bases]

def filter_reverse_complements(sequences):
    """
    From a list of sequences, keep only one of each reverse‑complement pair
    to avoid redundant scoring.
    """
    final = []
    seen = set()
    for seq in sequences:
        rc = reverse_complement(seq)
        if seq not in seen and rc not in seen:
            final.append(seq)
            seen.add(seq)
    return final

def refine_test_set_for_constraints(test_set, forbidden, fixed):
    """
    Remove any sequences (or their complements) that are in the forbidden
    or fixed lists, returning the allowable pool for sampling.
    """
    out = [s for s in test_set if s not in forbidden and s not in fixed]
    for seq in forbidden + fixed:
        rc = reverse_complement(seq)
        if rc in out:
            out.remove(rc)
    return out

def compute_target_scores(seq, combined, potapov):
    """
    For a given overhang `seq`, compute:
      - 'on'‑target score: interactions with its reverse complement
      - 'off'‑target score: interactions with all other candidates
    using the Potapov empirical dataset.
    """
    rc = reverse_complement(seq)
    row_seq = potapov.loc[potapov["Overhang"] == seq].squeeze()
    row_rc  = potapov.loc[potapov["Overhang"] == rc ].squeeze()
    on  = row_seq.get(rc, 0) + row_rc.get(seq, 0)
    others = [o for o in combined if o not in (seq, rc)]
    off = row_seq.loc[others].sum() + row_rc.loc[others].sum()
    return on, off

def calculate_fidelity_score(seqs, potapov):
    """
    Compute overall fidelity for a set of overhangs as the product over
    each seq of (1 – off/(on+off)), ensuring non‑redundant scoring.
    """
    seqs = filter_reverse_complements(seqs)
    combined = list(set(seqs + [reverse_complement(s) for s in seqs]))
    scores = []
    for s in seqs:
        on, off = compute_target_scores(s, combined, potapov)
        scores.append(1 - off/(on + off))
    return np.prod(scores)

def should_accept(current, new, temp):
    """
    Metropolis acceptance: always accept if `new` > `current`,
    otherwise accept with probability exp((new-current)/temp).
    """
    if new > current:
        return True
    return math.exp((new-current)/temp) > random.random()

def mutate(current, pool, fixed):
    """
    Propose a new candidate by swapping out one non‑fixed overhang
    for another from the allowed pool.
    """
    choices = [s for s in current if s not in fixed]
    to_replace = random.choice(choices)
    new = random.choice([s for s in pool if s not in current and s not in fixed])
    return [new if s == to_replace else s for s in current]

def monte_carlo_optimization(pool, set_size, iterations, start_temp, cooling,
                             fixed, conv_thres, potapov, initial=None):
    """
    Perform simulated annealing:
      1. Initialize candidate set (optionally seeded)
      2. Loop for `iterations`:
         - mutate, score, decide acceptance
         - track global best and a top‑50 list
         - cool temperature, early‑stop if plateaued
    Returns: best_set, best_score, list_of_top_sets.
    """
    # Initialize starting set
    if initial:
        current = initial
    else:
        available = [s for s in pool if s not in fixed]
        current = fixed + random.sample(available, set_size - len(fixed))
    current_score = calculate_fidelity_score(current, potapov)
    best, best_score = current, current_score
    temp = start_temp
    no_imp = 0
    top_sets = []

    for _ in range(iterations):
        cand = mutate(current, pool, fixed)
        score = calculate_fidelity_score(cand, potapov)

        # Accept or reject
        if should_accept(current_score, score, temp):
            current, current_score = cand, score

            # Update best if improved
            if score > best_score:
                improvement = (score - best_score) / best_score
                no_imp = 0 if improvement > conv_thres else no_imp + 1
                best, best_score = cand, score
            else:
                no_imp += 1

        # Cool down
        temp *= cooling

        # Early stop if stuck
        if no_imp >= 1500:
            print("Stopped: no improvement.")
            break

        # Maintain top‑50 list
        if len(top_sets) < 50:
            top_sets.append((cand, score))
        else:
            lowest = min(top_sets, key=lambda x: x[1])
            if score > lowest[1]:
                top_sets.remove(lowest)
                top_sets.append((cand, score))

    top_sets.sort(key=lambda x: x[1], reverse=True)
    return best, best_score, top_sets

# ----------------------------------------
# Main entry point
# ----------------------------------------
def main():
    """
    Parse arguments, load data, build candidate pool, optionally seed
    from hinge‑sets, run Monte Carlo optimizer, and save results.
    """
    parser = argparse.ArgumentParser(
        description='Monte Carlo Optimization for Sequence Selection.'
    )
    parser.add_argument('--potapov_data',       type=str,
                        default='data/FileS03_T4_18h_25C.xlsx')
    parser.add_argument('--set_size',           type=int,    default=30)
    parser.add_argument('--iterations',         type=int,    default=50000)
    parser.add_argument('--start_temperature',  type=float,  default=100)
    parser.add_argument('--cooling_rate',       type=float,  default=0.995)
    parser.add_argument('--fixed_overhangs',    type=str,
                        default='AATG,GCTT',
                        help='Comma‑separated fixed overhangs.')
    parser.add_argument('--forbidden_overhangs', type=str, nargs='*',
                        default=[
                            'GGAC','CGCG','GGCC','GCCA','CCGC','GGCG','GTCG','TCGC',
                            'GGGT','GGGG','TGCG','TAAA','GCAT','GGCT','TTTA','GCGG',
                            'TATA','GGAT','CCCC','GGGC','GCGT','TTAA','GGAG','GGTG',
                            'ATAT','CGCG','GCGC','GATC','CTAG','GTAC','CATG','AATT',
                            'TTAA','CCGG','GGCC','AGCT','TCGA','ACGT','TGCA'
                        ],
                        help='Space or comma‑separated forbidden overhangs.')
    parser.add_argument('--test_set_type',      type=str,
                        choices=['with_palindromes','without_palindromes'],
                        default='without_palindromes')
    parser.add_argument('--convergence_threshold', type=float,
                        default=0.0005)
    parser.add_argument('--output_file',        type=str,    default=None)
    parser.add_argument('--use_hingesets',      action='store_true')
    args = parser.parse_args()

    # Parse fixed_overhangs into a list
    fixed_list = [s.strip() for s in args.fixed_overhangs.split(',') if s.strip()]

    # Parse forbidden_overhangs (allow comma‑separated single string)
    raw = args.forbidden_overhangs
    if len(raw) == 1 and ',' in raw[0]:
        forb_list = [s.strip() for s in raw[0].split(',') if s.strip()]
    else:
        forb_list = raw

    print("Fixed:    ", fixed_list)
    print("Forbidden:", forb_list)

    # Load Potapov interaction data
    pot = pd.read_excel(args.potapov_data)

    # Build initial candidate pool of 4‑mers (optionally excluding palindromes)
    all4 = generate_all_4mers()
    uniq = filter_reverse_complements(all4)
    if args.test_set_type == 'with_palindromes':
        test = uniq
    else:
        pals = [s for s in all4 if s == reverse_complement(s)]
        test = [s for s in uniq if s not in pals]
    test = refine_test_set_for_constraints(test, forb_list, fixed_list)

    # Optionally seed with top sets from a precomputed hingesets.xlsx
    initial = None
    if args.use_hingesets:
        try:
            df_h = pd.read_excel("data/hingesets.xlsx")
            df_h['set_size'] = df_h['set_size'].astype(int)
            df_h = df_h[df_h['set_size'] == args.set_size]
            if not df_h.empty:
                elites = [
                    [s.strip() for s in row['hf_oh_set'].split(',')]
                    for _, row in df_h.sort_values('fidelity', ascending=False).head(5).iterrows()
                ]
                initial = random.choice(elites)
                print("Seeding from hingesets:", initial)
        except Exception as e:
            print("Hingesets load error:", e)

    # Run the Monte Carlo optimizer
    best_set, best_score, top_sets = monte_carlo_optimization(
        test,
        args.set_size,
        args.iterations,
        args.start_temperature,
        args.cooling_rate,
        fixed_list,
        args.convergence_threshold,
        pot,
        initial
    )

    # Prepare run parameters for output
    params = {
        'potapov_data':   args.potapov_data,
        'set_size':       args.set_size,
        'iterations':     args.iterations,
        'start_temp':     args.start_temperature,
        'cooling_rate':   args.cooling_rate,
        'fixed':          ','.join(fixed_list),
        'forbidden':      ','.join(forb_list),
        'conv_thres':     args.convergence_threshold,
        'test_set_type':  args.test_set_type,
        'use_hingesets':  args.use_hingesets
    }
    df_params = pd.DataFrame(list(params.items()),
                             columns=['Parameter','Value'])

    # Determine output filename and write Excel workbook
    ts = datetime.now().strftime('%Y%m%d_%H%M%S')
    base = args.output_file or f"mc_{args.set_size}_{args.iterations}_{ts}"
    out_file = f"{base}.xlsx"
    with pd.ExcelWriter(out_file, engine='openpyxl') as writer:
        pd.DataFrame(top_sets, columns=['Set','Score']) \
          .to_excel(writer, sheet_name='TopSets', index=False)
        df_params.to_excel(writer, sheet_name='Parameters', index=False)

    # Report results
    print(f"Best score: {best_score:.4f}")
    print(f"Best set: {best_set}")
    print(f"Results saved to {out_file}")

if __name__ == '__main__':
    main()
