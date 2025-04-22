#!/usr/bin/env python3

"""
This script implements a protocol for gene library synthesis using golden gate assembly, 
similar in principle to the OMEGA approach descrbed the Romero lab 
(https://doi.org/10.1101/2025.03.22.644747 & https://github.com/RomeroLab/omega).

A fasta file with target coding sequences is  split into batches of roughly equal 
total length; the sequences are purged of default IIS sites and optionally codon-optimized.
Each batch is broken into indexed oligos that can be reassembled after PCR using the 
index primer pair for that batch. A Monte Carlo optimization approach is employed to
help identfiy high-fidelity overhang sets for each batch.

Usage:
    python iggypop/multiplexed_run.py --i test/35_.fasta --o output_prefix --yml config.yml [options]
"""

import os
import io
import sys
import argparse
import yaml
import random
import itertools
import shutil
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from concurrent.futures import ProcessPoolExecutor, as_completed
import contextlib

# Core hinging logic: generate and evaluate overhang sets
from chisel_hinge import get_overhang_sets, find_cut_solution, chisel

# Helper functions for logging, fidelity calculations, and sequence metrics
from pop_helpers import (
    log_and_print,
    calculate_fidelity_score,
    calculate_segment_length,
    calculate_codon_frequencies,
    reverse_complement
)

# Fallback codon table loader
from python_codon_tables import get_codons_table


# -----------------------------------------------------------------------------
# Utility Functions
# -----------------------------------------------------------------------------

def parse_fasta_lengths(fasta_path):
    """
    Read a FASTA file and return:
      - lengths: list of (sequence_id, sequence_length)
      - sequences: dict mapping sequence_id to sequence string
    """
    lengths, sequences = [], {}
    for rec in SeqIO.parse(fasta_path, 'fasta'):
        lengths.append((rec.id, len(rec.seq)))
        sequences[rec.id] = str(rec.seq)
    return lengths, sequences


import math

def assign_to_batches(lengths, max_batch_size):
    """
    Balance gene‐counts across the minimum number of bins
    needed to respect max_batch_size.
    """
    # 1) Compute the minimum number of bins by total length
    total_len = sum(l for _, l in lengths)
    n_bins = max(1, math.ceil(total_len / max_batch_size))

    # 2) Initialize bins, tracking both sum and count
    batches = [[] for _ in range(n_bins)]
    sums    = [0] * n_bins
    counts  = [0] * n_bins

    # 3) Place longest genes first, into the bin with
    #    (a) smallest count, then (b) smallest current sum
    for seq_id, length in sorted(lengths, key=lambda x: x[1], reverse=True):
        # find best bin that can still fit this gene
        idx = min(
            (i for i in range(n_bins) if sums[i] + length <= max_batch_size),
            key=lambda i: (counts[i], sums[i]),
            default=None
        )
        if idx is None:
            # if none fit, open a new bin
            idx = n_bins
            batches.append([])
            sums.append(0)
            counts.append(0)
            n_bins += 1

        # assign
        batches[idx].append(seq_id)
        sums[idx]   += length
        counts[idx] += 1

    return batches, sums


def load_yaml(yaml_path, log_file, cli_args):
    """
    Load default parameters from YAML, override with any CLI args provided,
    and log the final configuration.

    Returns:
      - data: dict of final parameters
    """
    with open(yaml_path) as f:
        data = yaml.safe_load(f)

    # Map CLI arguments to YAML keys
    overrides = {
        'primer_index': cli_args.primer_index,
        'n_gene_orders': cli_args.n_gene_orders,
        'swap_cycles': cli_args.swap_cycles,
        'n_tries': cli_args.n_tries,
        'codon_tbl': cli_args.codon_tbl,
        'codon_opt': cli_args.codon_opt,
        'original_species': cli_args.original_species,
        'max_batch_size': cli_args.max_batch_size
    }

    # Apply overrides when CLI args are explicitly set
    for key, val in overrides.items():
        if val is not None:
            data[key] = val

    # Log the command and all parameters
    log_and_print(f"Command: {' '.join(sys.argv)}", log_file)
    for k, v in data.items():
        log_and_print(f"{k}: {v}", log_file, quiet="on")

    return data


# -----------------------------------------------------------------------------
# Core Hinging and Fragment Generation
# -----------------------------------------------------------------------------

def process_order(records, full_oh_sets, external_overhangs,
                  yaml_defaults, potapov_data, segment_length, log_file):
    """
    For a given ordering of gene records:
      - Prepend/appand base ends to each sequence
      - Run Monte Carlo search (find_cut_solution) to assign internal overhangs
      - Calculate per-gene and batch-level fidelity scores
      - Generate fragment entries (with left/right overhangs)

    Returns:
      - gene_results: list of dicts with per-gene fidelity and overhangs
      - fragments: list of dicts describing each fragment sequence
      - used_ohs: list of all internal overhangs chosen
      - overall_fid: total fidelity score for the batch
    """
    used_ohs, fragments, gene_results = [], [], []
    start_idx = 0
    total_sets = len(full_oh_sets)

    # Unpack cutting/adapter parameters from config
    base5 = yaml_defaults['base_5p_end']
    base3 = yaml_defaults['base_3p_end']
    p5cut = yaml_defaults['pcr_5p_cut']
    p3cut = yaml_defaults['pcr_3p_cut']
    base_tries = yaml_defaults['n_tries']

    for idx, rec in enumerate(records):
        # Add base adapters
        seq_full = base5 + str(rec.seq) + base3
        # Increase tries for later genes to explore more options
        n_tries = int(base_tries * (1.2 ** idx))

        # Search for the best cut solution in available overhang sets
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            best_sol, _, _, best_idx = find_cut_solution(
                sequence=seq_full,
                overhang_sets=full_oh_sets[start_idx:],
                radius=yaml_defaults['radius'],
                external_overhangs=external_overhangs + used_ohs,
                segment_length=segment_length,
                n_tries=n_tries,
                potapov_data=potapov_data
            )
        # If no solution found, backtrack start index and continue
        if not best_sol:
            start_idx = max(0, start_idx - 1)
            continue

        # Record chosen overhang sequences
        new_ohs = [d['sequence'] for d in best_sol]
        used_ohs.extend(new_ohs)

        # Compute fidelity for this gene
        fid = calculate_fidelity_score(new_ohs + external_overhangs, potapov_data)
        gene_results.append({
            'accession': rec.id,
            'individual_fidelity': fid,
            'internal_ohs': new_ohs
        })

        # Build fragment sequences based on cut locations
        cuts = sorted(d['location'] for d in best_sol)
        prev = 0
        for i, cut in enumerate(cuts):
            core = seq_full[prev:cut] + new_ohs[i]
            frag_seq = p5cut + core + p3cut
            fragments.append({
                'accession': rec.id,
                'frag_id': f"{rec.id}_FRAG_{i+1}",
                'sequence': frag_seq,
                'left_oh': new_ohs[i-1] if i > 0 else external_overhangs[0],
                'right_oh': new_ohs[i],
                'frag_length': len(frag_seq)
            })
            prev = cut

        # Handle final fragment after the last cut
        core = seq_full[prev:]
        frag_seq = p5cut + core + p3cut
        fragments.append({
            'accession': rec.id,
            'frag_id': f"{rec.id}_FRAG_{len(cuts)+1}",
            'sequence': frag_seq,
            'left_oh': new_ohs[-1],
            'right_oh': external_overhangs[-1],
            'frag_length': len(frag_seq)
        })

        # Update start index to avoid reusing the same OH sets immediately
        row_used = start_idx + best_idx
        if row_used + 1 >= total_sets:
            start_idx = max(0, total_sets - 100)
        else:
            start_idx = row_used + 1

    # Compute overall batch fidelity
    overall_fid = calculate_fidelity_score(used_ohs + external_overhangs, potapov_data)
    return gene_results, fragments, used_ohs, overall_fid


def run_pipeline(records, out_dir, yaml_defaults, potapov_data,
                 segment_length, full_oh_sets, external_overhangs,
                 index_df, index_start, n_gene_orders, swap_cycles, log_file):
    """
    Execute one batch:
      1) Global exploration over n_gene_orders random restarts
      2) Local exploitation via up to swap_cycles random swaps
      3) Attach index primers, write primer files, update fragments
      4) Return (gene_results, updated_fragments, used_overhangs, batch_fidelity, primers)
    """

    def attempt_trials(tries):
        best_res = None
        best_fid = -1
        best_order = None

        # 1) Global exploration: random restarts
        for t in range(tries):
            trial = records if t == 0 else random.sample(records, len(records))
            gr, fr, ohs, fid = process_order(
                trial,
                full_oh_sets,
                external_overhangs,
                yaml_defaults,
                potapov_data,
                segment_length,
                log_file
            )
            if fid > best_fid:
                best_fid, best_res, best_order = fid, (gr, fr, ohs), trial

        # 2) Local exploitation: random-swap hill-climbing
        if best_res:
            gr, fr, ohs = best_res
            n = len(best_order)
            max_swaps = min(swap_cycles, n * (n - 1) // 2)
            for a, b in random.sample(list(itertools.combinations(range(n), 2)), max_swaps):
                trial = best_order.copy()
                trial[a], trial[b] = trial[b], trial[a]
                gr2, fr2, ohs2, fid2 = process_order(
                    trial,
                    full_oh_sets,
                    external_overhangs,
                    yaml_defaults,
                    potapov_data,
                    segment_length,
                    log_file
                )
                if fid2 > best_fid:
                    best_fid, best_res, best_order = fid2, (gr2, fr2, ohs2), trial

        return best_res, best_fid

    # -------------------------------------------------------------------------
    # Run combined global + local search
    # -------------------------------------------------------------------------
    best_res, best_fid = attempt_trials(n_gene_orders)

    if not best_res:
        log_and_print(
            f"No successful solutions after {n_gene_orders} restarts + "
            f"hill‑climbing for batch {out_dir}.",
            log_file
        )
        return [], [], [], 0, []

    # Unpack best result
    (gr, fr, ohs), batch_fid = best_res, best_fid

    # -------------------------------------------------------------------------
    # Primer selection & writing
    # -------------------------------------------------------------------------
    row = index_df.iloc[index_start]
    primer_name = row['name']
    F_seq = row['F_seq']
    R_seq = row.get('R_seq', reverse_complement(row.get('R_seq_rc', '')))

    # Annotate gene results
    for entry in gr:
        entry['index_primer'] = primer_name

    primers = [{
        'primer_name': primer_name,
        'F_seq': F_seq,
        'R_seq': R_seq
    }]

    # Write primer CSV
    pd.DataFrame(primers).to_csv(
        os.path.join(out_dir, f"{yaml_defaults['run_name']}_index_primers_required.csv"),
        index=False
    )
    # Write primer FASTA
    with open(os.path.join(out_dir, f"{yaml_defaults['run_name']}_index_primers_required.fasta"), 'w') as pf:
        pf.write(f">{primer_name}_F\n{F_seq}\n")
        pf.write(f">{primer_name}_R\n{R_seq}\n")

    # -------------------------------------------------------------------------
    # Attach primers to fragments
    # -------------------------------------------------------------------------
    updated_fr = []
    for frag in fr:
        frag['left_oh']   = frag.get('left_oh') or external_overhangs[0]
        frag['right_oh']  = frag.get('right_oh') or external_overhangs[-1]
        frag['primer_name'] = primer_name
        frag['sequence']  = F_seq + frag['sequence'] + R_seq
        updated_fr.append(frag)

    # Write fragments FASTA
    with open(os.path.join(out_dir, 'fragments.fasta'), 'w') as wf:
        for frag in updated_fr:
            wf.write(f">{frag['frag_id']}\n{frag['sequence']}\n")

    return gr, updated_fr, ohs, batch_fid, primers


# ----------------------------------
# Worker for a single batch
# ----------------------------------
def batch_worker(params):
    import os

    (bidx, ids, chiselled_records, base_out, yaml_defaults,
     potapov_data, segment_length, full_oh_sets,
     external_overhangs, index_df, base_index, args, log_path) = params

    # ensure log directory exists
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    # make batch directory
    bdir = os.path.join(base_out, f"batch_{bidx}")
    os.makedirs(bdir, exist_ok=True)

    # open and run
    with open(log_path, 'w') as log_f:
        gr, fr, ohs, bfid, primers = run_pipeline(
            [r for r in chiselled_records if r.id in ids],
            bdir,
            yaml_defaults,
            potapov_data,
            segment_length,
            full_oh_sets,
            external_overhangs,
            index_df,
            base_index + (bidx - 1),
            args.n_gene_orders,
            args.swap_cycles,
            log_f
        )

    return bidx, gr, fr, ohs, bfid, primers


# ----------------------------------
# Main orchestration
# ----------------------------------
def main():
    # 1) Parse CLI arguments
    parser = argparse.ArgumentParser(
        description="Multiplexed hinging pipeline -- parallel batch execution"
    )
    parser.add_argument('--i', required=True, help="Input FASTA of CDS sequences")
    parser.add_argument('--o', required=True, help="Output prefix")
    parser.add_argument('--yml', default="yaml/domesticate_cds_minimal.yml", help="YAML config file")
    parser.add_argument('--primer_index', type=int, default=1,
                        help="1-based row in index CSV to start from")
    parser.add_argument('--max_batch_size', type=int, default=9000,
                        help="Max total bp per batch")
    parser.add_argument('--n_gene_orders', type=int, default=15,
                        help="Number of gene-order trials per batch")
    parser.add_argument('--swap_cycles', type=int, default=10,
                        help="Number of random swaps per batch")
    parser.add_argument('--n_tries', type=int, default=10,
                        help="Base number of hinging attempts per gene")
    parser.add_argument('--codon_opt', default="none",
                        help="Codon optimization method")
    parser.add_argument('--codon_tbl', default="cocoputs",
                        help="Codon table source")
    parser.add_argument('--original_species', help="Species for codon table")
    args = parser.parse_args()

    # 2) Setup output directories & log
    run_name = os.path.splitext(os.path.basename(args.o))[0]
    base_out = os.path.join('out', run_name)
    os.makedirs(base_out, exist_ok=True)
    log_file = open(os.path.join(base_out, 'run.log'), 'w')

    # Copy assets
    assets_dir = os.path.join(base_out, 'assets')
    os.makedirs(assets_dir, exist_ok=True)
    shutil.copy(args.i, assets_dir)
    shutil.copy(os.path.expanduser(args.yml), assets_dir)
    for fn in ['pop_helpers.py', 'chisel_hinge.py', 'multiplexed_run.py']:
        src = os.path.join('iggypop', fn)
        if os.path.exists(src):
            shutil.copy(src, assets_dir)

    # 3) Load & finalize config
    yaml_defaults = load_yaml(os.path.expanduser(args.yml), log_file, args)
    yaml_defaults.update({
        'n_tries': args.n_tries or yaml_defaults['n_tries'],
        'codon_tbl': args.codon_tbl or yaml_defaults.get('codon_tbl', 'cocoputs'),
        'codon_opt': args.codon_opt or yaml_defaults.get('codon_opt', 'none'),
        'original_species': args.original_species
                             or yaml_defaults.get('original_species')
                             or yaml_defaults.get('species', 'none'),
        'reports_dir': os.path.join(base_out, 'reports'),
        'run_name': run_name
    })
    os.makedirs(yaml_defaults['reports_dir'], exist_ok=True)

    # 4) Prepare codon table
    if yaml_defaults['codon_tbl'] == 'cocoputs':
        codon_tbl_obj, _ = calculate_codon_frequencies(
            'data/cleaned_coco.tsv',
            yaml_defaults['original_species']
        )
    else:
        codon_tbl_obj = get_codons_table(yaml_defaults['original_species'])
    yaml_defaults['codon_tbl_obj'] = codon_tbl_obj

    # 5) Load support data
    potapov_data   = pd.read_excel(yaml_defaults['fidelity_data'])
    segment_length = calculate_segment_length(
        yaml_defaults['pcr_5p_cut'],
        yaml_defaults['primer_length'],
        yaml_defaults['oligo_length']
    )
    external_overhangs = yaml_defaults['ext_overhangs']
    full_oh_sets      = get_overhang_sets(yaml_defaults['ohsets'], external_overhangs)
    index_df          = pd.read_csv(yaml_defaults['index_primers'])
    base_index        = args.primer_index - 1

    # 6) Chisel (codon-optimize)
    raw_records = list(SeqIO.parse(args.i, 'fasta'))
    chiselled_records = []
    for record in raw_records:
        print(f"Domesticating/optimizing {record.id} per {os.path.basename(args.yml)} specs")
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            optimized = chisel(
                str(record.seq),
                yaml_defaults['codon_opt'],
                codon_tbl_obj,
                yaml_defaults['original_species'],
                yaml_defaults.get('intron_constraints', []),
                yaml_defaults.get('loc_constraints_left', []),
                yaml_defaults.get('loc_constraints_right', []),
                yaml_defaults.get('left_bounds', []),
                yaml_defaults.get('right_bounds', []),
                log_file,
                True,
                yaml_defaults.get('pct', 0),
                0,
                0,
                quiet='on',
                file=os.path.expanduser(args.yml),
                results_path=os.path.join(
                    yaml_defaults['reports_dir'],
                    f"{record.id}.html"
                )
            )
        if optimized:
            record.seq = Seq(optimized)
            chiselled_records.append(record)
        else:
            log_and_print(f"Chiseling failed for {record.id}, skipping.", log_file)

    # Map of final sequences
    chiseled_map = {rec.id: str(rec.seq) for rec in chiselled_records}

    # Write designed full-length FASTA
    designed_fa = os.path.join(base_out, f"{run_name}_designed_seqs.fasta")
    with open(designed_fa, 'w') as df:
        for rec in chiselled_records:
            df.write(
                f">{rec.id}.chisel\n"
                f"{yaml_defaults['base_5p_end']}{rec.seq}{yaml_defaults['base_3p_end']}\n"
            )
    log_and_print(f"Wrote chiseled seqs to {designed_fa}", log_file)

    # 7) Split into batches (manual greedy)
    # get original lengths and sequences
    orig_lengths_list, orig_seqs = parse_fasta_lengths(args.i)
    orig_lens = dict(orig_lengths_list)
    # filter to chiselled
    id_lens = [(rec.id, orig_lens[rec.id]) for rec in chiselled_records]
    batches = []
    cur_batch, cur_len = [], 0
    for gid, glen in id_lens:
        if cur_batch and (cur_len + glen > args.max_batch_size):
            batches.append(cur_batch)
            cur_batch, cur_len = [], 0
        cur_batch.append(gid)
        cur_len += glen
    if cur_batch:
        batches.append(cur_batch)

    print("Partitioning into batches:")
    print("Batch\tGenes\tSize(bp)")
    for idx, batch in enumerate(batches, start=1):
        total_bp = sum(orig_lens[g] for g in batch)
        print(f"{idx}\t{len(batch)}\t{total_bp}")

    # 8) Parallel execution
    batch_params = []
    for bidx, batch in enumerate(batches, start=1):
        log_path = os.path.join(base_out, f"batch_{bidx}", "batch.log")
        batch_params.append((
            bidx, batch, chiselled_records, base_out, yaml_defaults,
            potapov_data, segment_length, full_oh_sets,
            external_overhangs, index_df, base_index, args, log_path
        ))

    all_primers, all_gene_data, all_frag_data, failed_batches = [], [], [], []
    with ProcessPoolExecutor() as executor:
        for bidx, gr, fr, ohs, bfid, primers in executor.map(batch_worker, batch_params):
            all_primers.extend(primers)
            if not gr:
                failed_batches.append(bidx)
                continue

            num_overhangs      = len(ohs)
            batch_total_length = sum(orig_lens[g['accession']] for g in gr)

            for g in gr:
                all_gene_data.append({
                    'accession':           g['accession'],
                    'pop_id':              f"{g['accession']}.chisel",
                    'i_seq':               orig_seqs[g['accession']],
                    'chiseled_seq':        f"{yaml_defaults['base_5p_end']}{chiseled_map[g['accession']]}{yaml_defaults['base_3p_end']}",
                    'index_primer':        g['index_primer'],
                    'genes_per_batch':     len(batch),
                    'batch_number':        bidx,
                    'individual_fidelity': g['individual_fidelity'],
                    'batch_fidelity':      bfid,
                    'num_overhangs':       num_overhangs,
                    'batch_total_length':  batch_total_length
                })

            for frag in fr:
                popid    = f"{frag['accession']}.1"
                suffix   = frag['frag_id'].split(frag['accession'])[1]
                left_oh  = frag.get('left_oh')  or external_overhangs[0]
                right_oh = frag.get('right_oh') or external_overhangs[-1]
                all_frag_data.append({
                    'accession':    frag['accession'],
                    'frag_id':      f"{popid}{suffix}",
                    'primer_index': frag['primer_name'],
                    'oligo':        frag['sequence'],
                    'left_oh':      left_oh,
                    'right_oh':     right_oh,
                    'frag_length':  frag['frag_length']
                })

    # 9) Final aggregation & outputs
    pool_fa = os.path.join(base_out, f"{run_name}_oligo_pool.fasta")
    with open(pool_fa, 'w') as pf:
        for e in all_frag_data:
            pf.write(f">{e['frag_id']}\n{e['oligo']}\n")
    log_and_print(f"Wrote oligo pool to {pool_fa}", log_file)

    primers_df = pd.DataFrame([
        {'name': p['primer_name'], 'F_seq': p['F_seq'], 'R_seq': p['R_seq']}
        for p in all_primers
    ])
    idx_csv = os.path.join(base_out, f"{run_name}_index_primers_required.csv")
    primers_df.to_csv(idx_csv, index=False)
    idx_fa  = os.path.join(base_out, f"{run_name}_index_primers_required.fasta")
    with open(idx_fa, 'w') as pf:
        for _, row in primers_df.iterrows():
            pf.write(f">{row['name']}_F\n{row['F_seq']}\n")
            pf.write(f">{row['name']}_R\n{row['R_seq']}\n")
    log_and_print(f"Wrote aggregated index primers to {idx_csv}", log_file)

    gene_df = pd.DataFrame(all_gene_data).drop(columns=['accession_fidelity'], errors='ignore')
    if not gene_df.empty:
        cols = ['accession','pop_id'] + [c for c in gene_df.columns if c not in ('accession','pop_id')]
        gene_df = gene_df[cols]
    frag_df = pd.DataFrame(all_frag_data)
    excel_path = os.path.join(base_out, f"{run_name}_all_data.xlsx")
    with pd.ExcelWriter(excel_path) as writer:
        frag_df.to_excel(writer,'fragment_data',index=False)
        if not gene_df.empty:
            gene_df.to_excel(writer,'seq_data',index=False)
    log_and_print(f"Wrote all data to {excel_path}", log_file)

    if failed_batches:
        log_and_print(
            f"WARNING: No valid solutions in batches {failed_batches}. Consider increasing --n_gene_orders or removing problematic genes.",
            log_file
        )
    log_file.close()


if __name__ == '__main__':
    main()
