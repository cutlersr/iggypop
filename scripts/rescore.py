#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import os
import re
from concurrent.futures import ProcessPoolExecutor

# Purpose:
#  - Calculate overhang set fidelity using specified scoring datasets.
#  - Compute new fidelities for T4_18h_37C, T4_BsmBI, T4_BsaI.
#  - Rename and reorder columns, add external_overhangs.
#  - Always write full rescored to top_rescored.xlsx and thinned subset to hingests.xlsx.
#  - Thinning uses farthest‑point sampling per set_size, default k=50.
#  - Remove exact duplicate sets before thinning and ensure each set appears only once.

SCORE_FILES = [
    'FileS04_T4_18h_37C.xlsx',
    'BsmBI-HFv2_T4.xlsx',
    'BsaI-HFv2_T4.xlsx'
]
FID_NAMES = ['fid_T4_18h_37c', 'fid_T4_BsmBI', 'fid_T4_BsaI']

def reverse_complement(seq: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(b, b) for b in reversed(seq))

def compute_target_scores(seq, combined, df):
    rc = reverse_complement(seq)
    row_seq = df.loc[df['Overhang'] == seq].squeeze()
    row_rc  = df.loc[df['Overhang'] == rc].squeeze()
    on = row_rc.get(seq, 0) + row_seq.get(rc, 0)
    valid_cols = set(df.columns) - {'Overhang'}
    others = [o for o in combined if o not in (seq, rc) and o in valid_cols]
    off = row_rc[others].sum() + row_seq[others].sum()
    return on, off

def fidelity_for_set(seqs, potapov_dfs):
    clean = [s.strip().upper() for s in seqs if all(c in 'ATCG' for c in s.strip().upper())]
    combined = list(set(clean + [reverse_complement(s) for s in clean]))
    scores = []
    for df in potapov_dfs:
        indiv = []
        for s in clean:
            on, off = compute_target_scores(s, combined, df)
            indiv.append(1 - off / (on + off))
        scores.append(np.prod(indiv))
    return scores

def parallel_scores(sets, dfs):
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(fidelity_for_set, s, dfs) for s in sets]
        return [f.result() for f in futures]

# ---- thinning utilities ----

def jaccard_distance(s1: set, s2: set) -> float:
    inter = s1 & s2
    union = s1 | s2
    return 1 - len(inter) / len(union) if union else 0

def farthest_point_sampling(sets: list, k: int) -> list:
    n = len(sets)
    # compute pairwise distances
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            d = jaccard_distance(sets[i], sets[j])
            dist[i, j] = dist[j, i] = d
    # seed with index 0 (highest fidelity after sorting)
    selected = [0]
    while len(selected) < min(k, n):
        best_idx, best_dist = None, -1
        for i in range(n):
            if i in selected: 
                continue
            d_min = min(dist[i, j] for j in selected)
            if d_min > best_dist:
                best_dist, best_idx = d_min, i
        selected.append(best_idx)
    return selected

# ---- main workflow ----

def main(input_file, data_folder, field,
         extra_files, do_thin, k, external_overhangs):

    # --- validations ---
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file '{input_file}' not found.")
    if not os.path.isdir(data_folder):
        raise NotADirectoryError(f"Data folder '{data_folder}' not found.")

    # --- read and check ---
    df = pd.read_excel(input_file)
    if 'raw_fidelity' not in df.columns:
        raise KeyError("Column 'raw_fidelity' missing; cannot proceed.")

    # duplicate raw fidelity
    df['fid_T4_18h_25c'] = df['raw_fidelity']

    # split sets
    if field not in df.columns:
        raise KeyError(f"Column '{field}' not in input; use --field to specify correct column.")
    sets = df[field].dropna().str.split(r',\s*')

    # load scoring datasets
    potapov_dfs, loaded = [], []
    for fname in SCORE_FILES + extra_files:
        path = os.path.join(data_folder, fname)
        if not os.path.exists(path):
            print(f"Warning: scoring file '{fname}' not found; skipping.")
            continue
        potapov_dfs.append(pd.read_excel(path))
        loaded.append(fname)
    print(f"Loaded scoring files: {loaded}")

    if len(potapov_dfs) < len(FID_NAMES):
        raise RuntimeError(f"Expected at least {len(FID_NAMES)} scoring datasets, got {len(potapov_dfs)}.")

    # compute new fidelity scores
    all_scores = parallel_scores(sets, potapov_dfs)
    for idx, name in enumerate(FID_NAMES):
        df[name] = [
            sl[idx] if idx < len(sl) else np.nan
            for sl in all_scores
        ]
    for col in FID_NAMES:
        df[col] = (df[col] * 100).round(1)

    # sort by primary fidelity
    df = df.sort_values(by='fid_T4_18h_25c', ascending=False)

    # rename columns
    df = df.rename(columns={
        'raw_fidelity': 'fidelity',
        field:          'hf_oh_set'
    })

    # round the fidelity column
    df['fidelity'] = df['fidelity'].round(1)

    # add external_overhangs
    df['external_overhangs'] = external_overhangs

    # drop unwanted columns
    df = df.drop(columns=['file', 'Score'], errors=True)

    # reorder columns
    final_cols = [
        'set_size',
        'hf_oh_set',
        'external_overhangs',
        'fidelity',
        'fid_T4_18h_25c',
        'fid_T4_18h_37c',
        'fid_T4_BsmBI',
        'fid_T4_BsaI'
    ]
    df = df[final_cols]

    # remove exact duplicates by set_size + hf_oh_set
    df = df.drop_duplicates(subset=['set_size', 'hf_oh_set'])

    # --- save full rescored ---
    full_out = 'top_rescored.xlsx'
    df.to_excel(full_out, index=False)
    print(f"Full rescored data saved to {full_out}")

    # --- thinning ---
    if do_thin:
        # parse hf_oh_set into Python sets once
        df['__oh_set__'] = df['hf_oh_set'].apply(lambda s: set(re.split(r',\s*', s)))
        pieces = []
        for size, group in df.groupby('set_size'):
            grp = group.sort_values('fidelity', ascending=False).reset_index(drop=True)
            # always sample up to k (min(k, n))
            idxs = farthest_point_sampling(list(grp['__oh_set__']), k)
            pieces.append(grp.iloc[idxs])
        df_thin = pd.concat(pieces, ignore_index=True)
        # ensure global uniqueness of hf_oh_set, keep highest fidelity
        df_thin = df_thin.sort_values('fidelity', ascending=False) \
                         .drop_duplicates(subset=['hf_oh_set']) \
                         .drop(columns='__oh_set__')
    else:
        df_thin = df.copy()

    # --- save thinned ---
    thin_out = 'hingests.xlsx'
    df_thin.to_excel(thin_out, index=False)
    print(f"{'Thinned' if do_thin else 'Unthinned'} data saved to {thin_out}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Rescore overhang sets and write both full and thinned outputs"
    )
    parser.add_argument('--i',  required=True, help='Input Excel file')
    parser.add_argument('--d',  default='../../data', help='Data folder with scoring files')
    parser.add_argument('--field', default='setA', help='Column containing overhang sets')
    parser.add_argument('--additional_files', nargs='*', default=[], help='Extra scoring files')
    parser.add_argument('--no-thin', dest='do_thin', action='store_false',
                        help='Skip thinning; hingests.xlsx will be identical to top_rescored.xlsx')
    parser.add_argument('--k', type=int, default=50,
                        help='Number of representatives per set_size when thinning')
    parser.add_argument('--external', default='AATG, GCTT',
                        help='Comma‑separated external overhangs string')
    args = parser.parse_args()

    main(args.i, args.d, args.field,
         args.additional_files, args.do_thin, args.k, args.external)
