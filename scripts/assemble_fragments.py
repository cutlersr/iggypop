#!/usr/bin/env python3
import re
import argparse
from collections import defaultdict, Counter

# ---------------- FASTA helpers ----------------
def parse_fasta(path):
    """Yield (header, sequence) tuples; header excludes leading '>'."""
    with open(path, "r") as fh:
        header = None
        seq_chunks = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)

def count_headers(path):
    with open(path, "r") as fh:
        return sum(1 for line in fh if line.startswith(">"))

# ---------------- Header parsing ----------------
# Supported header variants:
#  1) legacy FRAG:   ACC.N_FRAG_<frag>_PRI_<pri>
#  2) generic FRAG:  <BASE>_FRAG_<frag>_PRI_<pri>
#  3) two-index FRAG:<BASE>_FRAG_<block>_<frag>_PRI_<pri>  -> groups by (BASE, block)
#  4) L0:            <BASE>_L0_<block>_<frag>_PRI_<pri>
#  5) FINAL-ish:     "<name> ind_####"
#  6) cph tiles:     cph_fragm(e)?nts_ind_<NNNN>.<frag>    -> groups by stem before dot
#  7) pipe form:     "<anything>|index_primer_name=ind_####" -> treat as FINAL, name = left side
FRAG_2IDX_RE = re.compile(r"^(?P<base>.+?)_FRAG_(?P<blk>\d+)_(?P<frag>\d+)_PRI_(?P<pri>.+)$")
FRAG_1IDX_DOTTED_RE = re.compile(r"^(?P<acc>.+?)\.(?P<n>\d+)_FRAG_(?P<frag>\d+)_PRI_(?P<pri>.+)$")
FRAG_1IDX_GENERIC_RE = re.compile(r"^(?P<base>.+?)_FRAG_(?P<frag>\d+)_PRI_(?P<pri>.+)$")
L0_RE   = re.compile(r"^(?P<acc>.+?)_L0_(?P<block>\d+)_(?P<frag>\d+)_PRI_(?P<pri>.+)$")
FINAL_RE= re.compile(r"^(?P<name>.+?)\s+ind_\d+$", re.IGNORECASE)
CPH_RE  = re.compile(r"^(?P<stem>cph_fragm(?:e)?nts_ind_\d+)\.(?P<frag>\d+)$", re.IGNORECASE)
PIPE_FINAL_RE = re.compile(r"^(?P<left>.+?)\|index_primer_name=ind_\d+$", re.IGNORECASE)

def classify_header(h):
    """
    Return (type, base_key, order_key, pri) where:
      type in {"FRAG","L0","FINAL","LIBRARY","UNKNOWN"}
      base_key groups fragments that assemble together
      order_key determines assembly order (int or tuple)
      pri is the PRI_* suffix if present
    """
    hl = h.lower()
    if hl.startswith("linker_") or "library" in hl:
        return ("LIBRARY", h, None, None)

    m = PIPE_FINAL_RE.match(h)
    if m:
        return ("FINAL", m.group("left"), None, None)

    m = CPH_RE.match(h)
    if m:
        base = m.group("stem")
        order = int(m.group("frag"))
        return ("FRAG", base, order, None)

    # two-index FRAG: group BY BLOCK
    m = FRAG_2IDX_RE.match(h)
    if m:
        base = f"{m.group('base')}::blk{m.group('blk')}"   # separate groups per block
        order = int(m.group("frag"))
        return ("FRAG", base, order, m.group("pri"))

    # dotted legacy FRAG
    m = FRAG_1IDX_DOTTED_RE.match(h)
    if m:
        base = f"{m.group('acc')}.{m.group('n')}"
        order = int(m.group("frag"))
        return ("FRAG", base, order, m.group("pri"))

    # generic single-index FRAG
    m = FRAG_1IDX_GENERIC_RE.match(h)
    if m:
        base = m.group("base")
        order = int(m.group("frag"))
        return ("FRAG", base, order, m.group("pri"))

    # L0 two-step
    m = L0_RE.match(h)
    if m:
        base = f"{m.group('acc')}_L0_{m.group('block')}"
        order = int(m.group("frag"))
        return ("L0", base, order, m.group("pri"))

    # Already-assembled constructs with index primer suffix
    m = FINAL_RE.match(h)
    if m:
        name = m.group("name")
        return ("FINAL", name, None, None)

    return ("UNKNOWN", None, None, None)

# ---------------- Processing ----------------
def process(input_path, n_trim=25, skip_library=True):
    input_total = count_headers(input_path)

    ids = []
    recs = []
    for h, s in parse_fasta(input_path):
        ids.append(h)
        recs.append((h, s))

    id_counts = Counter(ids)
    dup_occ = sum(v-1 for v in id_counts.values())
    dup_ids = [k for k,v in id_counts.items() if v > 1]

    frag_groups = defaultdict(list)  # base -> list[(order, trimmed_seq)]
    l0_groups   = defaultdict(list)  # base -> list[(order, trimmed_seq)]
    finals      = {}                 # name -> seq
    unknown_ids = []

    lib_count = 0
    lib_len_sum = 0
    gene_count = 0
    gene_len_sum = 0
    full_len_sum = 0

    for h, s in recs:
        full_len_sum += len(s)
        typ, base, order, pri = classify_header(h)

        if typ == "LIBRARY":
            lib_count += 1
            trim = 26 if "_L0_" in h else n_trim
            ts = s[trim:-trim] if len(s) >= 2*trim else s
            lib_len_sum += len(ts)
            if not skip_library:
                pass
            continue

        if typ == "FRAG":
            trim = n_trim
            ts = s[trim:-trim] if len(s) >= 2*trim else s
            frag_groups[base].append((order, ts))
            gene_count += 1
            gene_len_sum += len(ts)
            continue

        if typ == "L0":
            trim = 26
            ts = s[trim:-trim] if len(s) >= 2*trim else s
            l0_groups[base].append((order, ts))
            gene_count += 1
            gene_len_sum += len(ts)
            continue

        if typ == "FINAL":
            finals[base] = s
            continue

        unknown_ids.append(h)

    return {
        "input_total": input_total,
        "records_total": len(recs),
        "unique_ids": len(id_counts),
        "dup_ids": dup_ids,
        "dup_occurrences": dup_occ,
        "frag_groups": frag_groups,
        "l0_groups": l0_groups,
        "finals": finals,
        "unknown_ids": unknown_ids,
        "gene_count": gene_count,
        "lib_count": lib_count,
        "gene_len_sum": gene_len_sum,
        "lib_len_sum": lib_len_sum,
        "full_len_sum": full_len_sum,
    }

# ---------------- Assembly ----------------
def assemble_groups(groups, cloning_5, cloning_3):
    """
    Assemble each base group by 4-bp overlaps.
    Emits separate contigs when a "reset" boundary is detected: prev[-4:]==cloning_3 and curr[:4]==cloning_5.
    Returns (assemblies, overhangs, mismatches), where
      assemblies: dict of id->sequence (id may include '::partN')
      overhangs:  id->list[4-mers] [5', internal..., 3']
      mismatches: id->list[str]
    """
    assemblies = {}
    overhangs  = {}
    mismatches = defaultdict(list)

    for base, items in groups.items():
        items_sorted = sorted(items, key=lambda x: (x[0] is None, x[0]))
        seqs = [s for _, s in items_sorted]
        if not seqs:
            continue

        part_idx = 1
        assembled = seqs[0]
        oh = [assembled[:4]]
        out_id = base if part_idx == 1 else f"{base}::part{part_idx}"

        for i in range(1, len(seqs)):
            prev = seqs[i-1]
            curr = seqs[i]
            prev_tail = prev[-4:]
            curr_head = curr[:4]
            # reset boundary? (e.g., GCTT -> AATG or TAGT/TCGC, depending on args)
            if prev_tail == cloning_3 and curr_head == cloning_5:
                # flush current contig
                oh.append(assembled[-4:])
                assemblies[out_id] = assembled
                overhangs[out_id] = oh
                part_idx += 1
                out_id = base if part_idx == 1 else f"{base}::part{part_idx}"
                # start new contig
                assembled = curr
                oh = [curr_head]
                continue

            # normal join
            oh.append(curr_head)
            if prev_tail != curr_head:
                mismatches[out_id].append(f"Mismatch {i-1}->{i}: {prev_tail} != {curr_head}")
            assembled = assembled[:-4] + curr

        # flush last contig
        oh.append(assembled[-4:])
        assemblies[out_id] = assembled
        overhangs[out_id] = oh

    return assemblies, overhangs, mismatches

# ---------------- Output ----------------
def write_fasta(records, output_file):
    with open(output_file, "w") as out:
        for rec_id, seq in records:
            out.write(f">{rec_id}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")

# ---------------- Main ----------------
def main():
    ap = argparse.ArgumentParser(description="Iggypop fragment assembler & checker (v6)")
    ap.add_argument("--i", required=True, help="Input FASTA")
    ap.add_argument("--o", default="assembled_seq.fasta", help="Output FASTA")
    ap.add_argument("--n", type=int, default=25, help="Trim N for FRAG_* (default 25). L0_* uses 26.")
    ap.add_argument("--cloning_ohs", nargs=2, default=["AATG", "GCTT"], help="Expected cloning overhangs 5'/3'")
    ap.add_argument("--skip_library", action="store_true", default=True, help="Skip assembling library/linker sequences")
    args = ap.parse_args()

    state = process(args.i, n_trim=args.n, skip_library=args.skip_library)

    # Assemble FRAG and L0 step1 (with reset-aware stitching)
    fr_assemblies, fr_overhangs, fr_mismatch = assemble_groups(state["frag_groups"], args.cloning_ohs[0], args.cloning_ohs[1])
    l0_step1_assemblies, l0_step1_overhangs, l0_mismatch = assemble_groups(state["l0_groups"], args.cloning_ohs[0], args.cloning_ohs[1])

    # Step 2: join L0 blocks per accession (group blocks by prefix before _L0_)
    groups_for_step2 = defaultdict(list)  # gene -> list of (block_num, seq)
    for base, seq in l0_step1_assemblies.items():
        m = re.search(r"_L0_(\d+)", base)
        block_num = int(m.group(1)) if m else 0
        gene = re.sub(r"_L0_.*", "", base)
        groups_for_step2[gene].append((block_num, seq))

    l0_final = {}
    l0_final_overhangs = {}
    step2_mismatch = []

    for gene, lst in groups_for_step2.items():
        lst.sort(key=lambda x: x[0])
        if not lst:
            continue
        parts = []
        for _, seq in lst:
            core = seq[11:-11] if len(seq) >= 22 else seq  # typical BsmBI cores
            parts.append({"core": core, "five": core[:4], "three": core[-4:]})
        if not parts:
            continue
        assembled = parts[0]["core"]
        oh = [parts[0]["five"]]
        for i in range(1, len(parts)):
            prev = parts[i-1]
            curr = parts[i]
            oh.append(curr["five"])
            if prev["three"] != curr["five"]:
                step2_mismatch.append(f"{gene} block{i-1}->{i}: {prev['three']} != {curr['five']}")
            assembled = assembled + curr["core"][4:]
        # enforce cloning_ohs at ends for step2 outputs
        if not assembled.startswith(args.cloning_ohs[0]):
            assembled = args.cloning_ohs[0] + assembled
        if not assembled.endswith(args.cloning_ohs[1]):
            assembled = assembled + args.cloning_ohs[1]
        oh.append(parts[-1]["three"])
        l0_final[gene] = assembled
        l0_final_overhangs[gene] = oh

    # Build final outputs
    final_records = []
    final_overhangs = {}

    for k, seq in fr_assemblies.items():
        final_records.append((k, seq))
        final_overhangs[k] = fr_overhangs.get(k, [])

    for k, seq in l0_step1_assemblies.items():
        rid = f"STEP1:{k}"
        final_records.append((rid, seq))
        final_overhangs[rid] = l0_step1_overhangs.get(k, [])

    for k, seq in l0_final.items():
        rid = f"STEP2:{k}"
        final_records.append((rid, seq))
        final_overhangs[rid] = l0_final_overhangs.get(k, [])

    for k, seq in state["finals"].items():
        rid = f"FINAL:{k}"
        final_records.append((rid, seq))
        final_overhangs[rid] = [seq[:4], seq[-4:]]

    write_fasta(final_records, args.o)

    # ---------------- Reporting ----------------
    print("\n--- Assembly Statistics ---")
    print(f"Total oligos in input file (headers): {state['input_total']}")
    print(f"Total records parsed: {state['records_total']}")
    print(f"Unique header IDs: {state['unique_ids']}")
    if state["dup_occurrences"]:
        print(f"Duplicate headers (occurrences beyond first): {state['dup_occurrences']} (unique dup IDs: {len(state['dup_ids'])})")

    print(f"Total oligos for genes (FRAG + L0): {state['gene_count']}")
    print(f"Total oligos for libraries/linkers: {state['lib_count']}")
    avg_trimmed = (state['gene_len_sum'] + state['lib_len_sum']) / max(1, state['records_total'])
    avg_full    = state['full_len_sum'] / max(1, state['records_total'])
    print(f"Average bp per oligo (trimmed est.): {avg_trimmed:.0f} bp")
    print(f"Average bp per oligo (including primers): {avg_full:.0f} bp")

    step1_final_count = len([rid for rid,_ in final_records if rid.startswith('STEP1:')])
    step2_final_count = len([rid for rid,_ in final_records if rid.startswith('STEP2:')])
    onestep_count     = len([rid for rid,_ in final_records if not rid.startswith(('STEP1:','STEP2:','FINAL:'))])
    finals_count      = len([rid for rid,_ in final_records if rid.startswith('FINAL:')])

    total_final = len(final_records)
    total_bases = sum(len(seq) for _, seq in final_records)
    avg_final_length = total_bases / max(1, total_final)

    if step2_final_count > 0:
        bases_step1 = sum(len(seq) for rid, seq in final_records if rid.startswith("STEP1:"))
        bases_step2 = sum(len(seq) for rid, seq in final_records if rid.startswith("STEP2:"))
        bases_onestep = sum(len(seq) for rid, seq in final_records if not rid.startswith(("STEP1:","STEP2:","FINAL:")))
        print(f"Total step 1 (L0) genes assembled final: {step1_final_count}")
        print(f"Total step 2 genes assembled final: {step2_final_count}")
        print(f"Total one-step genes assembled final: {onestep_count}")
        print(f"Total number of bases synthesized (step 1): {bases_step1}")
        print(f"Total number of bases synthesized (step 2): {bases_step2}")
        print(f"Total number of bases synthesized (one-step): {bases_onestep}")
    else:
        print(f"Total assembled records (all): {total_final}")
        print(f"Total number of bases synthesized: {total_bases}")

    print(f"Average assembled sequence length: {avg_final_length:.0f} bp")

    # Summaries
    if state["unknown_ids"]:
        print(f"Warning: {len(state['unknown_ids'])} sequences did not match expected formats (FRAG/L0/FINAL/LIBRARY). Showing first 10:")
        for u in state["unknown_ids"][:10]:
            print(f"  - {u}")

    total_fr_mis = sum(len(v) for v in fr_mismatch.values())
    if total_fr_mis:
        print(f"Step1 FRAG overlap mismatches: {total_fr_mis} total. Showing first 10:")
        shown = 0
        for msgs in fr_mismatch.values():
            for m in msgs:
                print(f"  - {m}")
                shown += 1
                if shown >= 10: break
            if shown >= 10: break

    total_l0_mis = sum(len(v) for v in l0_mismatch.values())
    if total_l0_mis:
        print(f"L0 step1 overlap mismatches: {total_l0_mis} total. Showing first 10:")
        shown = 0
        for msgs in l0_mismatch.values():
            for m in msgs:
                print(f"  - {m}")
                shown += 1
                if shown >= 10: break
            if shown >= 10: break

    if step2_final_count and step2_mismatch:
        print(f"Step2 overlap mismatches: {len(step2_mismatch)} total. Showing first 10:")
        for m in step2_mismatch[:10]:
            print(f"  - {m}")

    # Cloning overhang audit only on STEP2 and FINAL outputs (skip FRAG one-step unless requested)
    cloning_oh_mismatches = []
    for rid, seq in final_records:
        if rid.startswith("STEP2:") or rid.startswith("FINAL:"):
            if not (seq.startswith(args.cloning_ohs[0]) and seq.endswith(args.cloning_ohs[1])):
                five = seq[:4]
                three = seq[-4:]
                ohs = final_overhangs.get(rid, [five, three])
                cloning_oh_mismatches.append((rid, ohs, five, three))
    if cloning_oh_mismatches:
        print(f"Cloning overhang mismatches: {len(cloning_oh_mismatches)} total. Showing first 10:")
        for rid, ohs, five, three in cloning_oh_mismatches[:10]:
            print(f"  - {rid} ohs:{ohs} (5'={five}, 3'={three})")

if __name__ == "__main__":
    main()
