#!/usr/bin/env python3
"""
Demux Saunders 2025 Perturb-Multi sgRNA capture FASTQs to per-cell sgRNA counts.

Chemistry (reverse-engineered from raw FASTQ, see saunders2025_fastq_reprocessing.md):

  R1 (101 nt):
    [0:16]  16-nt 10x cell barcode      (Flex 737K-fixed-rna-profiling whitelist)
    [16:28] 12-nt UMI
    [28:58] ~30 nt of SplintR_barcode 3' end + linker  (TACGTGCT linker around pos 40-50)
    [58:78] 20-nt PROTOSPACER            (DISCRIMINATING — match against mmc5)
    [78:..] sgRNA scaffold start (GTTT[CA]AGAG...)

  R2 (101 nt):
    [0:8]   8-nt RNA Probe Barcode     (BC001=ACTTTAGG, BC002=AACGGGAA,
                                          BC003=AGTAGGCT, BC004=ATGTTGAC — 10x Flex)
    [8:..]  0-3 phasing N + constant region GCTATGCTGTTTCCAGCTTAGCTCT
    [..-..] SplintR_barcode 5' end portion + 10x CS1

Output per SRR:
   <out_dir>/<srr>.counts.tsv.gz   columns: cell_barcode, sample_bc, sgrna_id, gene, n_reads, n_umis
   <out_dir>/<srr>.stats.json      summary stats
"""
import argparse, csv, gzip, json, os, sys, time
from collections import defaultdict, Counter

# Cell barcode whitelist path
WHITELIST_PATH = "/nfs/sw/easybuild/software/CellRanger/8.0.1/lib/python/cellranger/barcodes/737K-fixed-rna-profiling.txt.gz"
PROBE_BC = {"ACTTTAGG": "BC001", "AACGGGAA": "BC002", "AGTAGGCT": "BC003", "ATGTTGAC": "BC004"}

# Empirical positions in R1
PROTOSPACER_POS_CANDIDATES = (58, 57, 59, 56, 60)
PROTOSPACER_LEN = 20
CB_END = 16
UMI_END = 28
PROBE_BC_END = 8

BASES = ("A", "C", "G", "T")


def load_whitelist():
    s = set()
    with gzip.open(WHITELIST_PATH, "rt") as fh:
        for line in fh:
            s.add(line.strip())
    return s


def load_protospacers(mmc5_path):
    """Returns: ps_to_sgid (dict 20-mer → (sgrna_id, gene, CR))."""
    ps_to_sg = {}
    sg_meta = {}
    with open(mmc5_path) as fh:
        rd = csv.DictReader(fh)
        for i, row in enumerate(rd):
            ps = row["protospacer"].upper()
            if ps.startswith("G") and len(ps) == 21:
                ps = ps[1:]
            if len(ps) != PROTOSPACER_LEN:
                continue
            sgid = f"sg{i:04d}_{row['gene']}_{row['CR']}"
            if ps in ps_to_sg:
                # Two sgRNAs with the same protospacer — keep both (rare)
                continue
            ps_to_sg[ps] = sgid
            sg_meta[sgid] = (row["gene"], row["CR"], row["SplintR_barcode"])
    return ps_to_sg, sg_meta


def build_1mm_index(seq_set):
    """Hamming-1 expansion. Maps each variant → original (only if uniquely correctable)."""
    out = {}
    ambig = set()
    for s in seq_set:
        out[s] = s   # self
    for s in seq_set:
        for i in range(len(s)):
            for b in BASES:
                if b == s[i]:
                    continue
                v = s[:i] + b + s[i+1:]
                if v in seq_set:
                    continue
                if v in out:
                    ambig.add(v)
                else:
                    out[v] = s
    for v in ambig:
        del out[v]
    return out


def correct_cb(cb, wl, wl_1mm):
    """Returns the corrected cell barcode or None."""
    if cb in wl:
        return cb
    return wl_1mm.get(cb)


def correct_probe_bc(pbc):
    """Returns BC001..BC004 or None."""
    if pbc in PROBE_BC:
        return PROBE_BC[pbc]
    # Try Hamming-1
    for raw, name in PROBE_BC.items():
        ms = sum(a != b for a, b in zip(pbc, raw))
        if ms <= 1:
            return name
    return None


def find_protospacer(r1, ps_to_sg, ps_1mm):
    """Try canonical position 58 first; allow 1mm; fall back to nearby positions.

    ps_to_sg: dict[20-mer protospacer] -> sgid
    ps_1mm  : dict[1-mismatch variant]  -> corrected 20-mer protospacer
    """
    for pos in PROTOSPACER_POS_CANDIDATES:
        window = r1[pos:pos + PROTOSPACER_LEN]
        if len(window) != PROTOSPACER_LEN:
            continue
        if window in ps_to_sg:
            return ps_to_sg[window], pos
        corrected = ps_1mm.get(window)
        if corrected is not None and corrected in ps_to_sg:
            return ps_to_sg[corrected], pos
    return None, None


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--srr", required=True)
    ap.add_argument("--mmc5", default="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/competitor_analysis/Saunders/mmc5.csv")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--max-reads", type=int, default=0, help="0 = no limit; >0 = stop early for testing")
    ap.add_argument("--no-cb-1mm", action="store_true", help="Skip CB Hamming-1 correction (faster)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    t0 = time.time()

    print(f"[{time.strftime('%T')}] Loading whitelist + protospacers...", flush=True)
    wl = load_whitelist()
    print(f"  whitelist: {len(wl):,}")
    ps_to_sg, sg_meta = load_protospacers(args.mmc5)
    print(f"  protospacers: {len(ps_to_sg):,}")
    print(f"  sgRNAs total: {len(sg_meta):,}")

    print(f"[{time.strftime('%T')}] Building 1-mismatch indices...", flush=True)
    ps_1mm = build_1mm_index(set(ps_to_sg.keys()))
    print(f"  protospacer 1mm index: {len(ps_1mm):,}")

    if args.no_cb_1mm:
        wl_1mm = {b: b for b in wl}
    else:
        # CB 1mm index is too large in memory (35M entries × 16 chars), skip and just use exact
        wl_1mm = {b: b for b in wl}

    # Map protospacer-id -> (gene, CR) for output
    ps_to_meta = {sgid: sg_meta[sgid] for sgid in sg_meta}

    # Aggregate: (cell_barcode_16, sample_bc, sgrna_id) -> set(UMI), n_reads
    counts = defaultdict(lambda: [0, set()])

    n_total = 0
    n_cb_valid = 0
    n_probe_valid = 0
    n_proto_valid = 0
    n_all_valid = 0
    n_proto_pos = Counter()

    print(f"[{time.strftime('%T')}] Streaming {args.r1} + {args.r2}...", flush=True)

    fh1 = open_maybe_gz(args.r1)
    fh2 = open_maybe_gz(args.r2)

    while True:
        # Block of 4 lines per record
        h1 = fh1.readline()
        if not h1:
            break
        seq1 = fh1.readline().rstrip()
        plus1 = fh1.readline()
        qual1 = fh1.readline()
        h2 = fh2.readline()
        seq2 = fh2.readline().rstrip()
        plus2 = fh2.readline()
        qual2 = fh2.readline()

        n_total += 1
        if args.max_reads and n_total > args.max_reads:
            break

        # 1) Cell barcode
        cb = seq1[:CB_END]
        cb_corr = correct_cb(cb, wl, wl_1mm)
        if cb_corr is None:
            continue
        n_cb_valid += 1

        # 2) UMI
        umi = seq1[CB_END:UMI_END]
        if "N" in umi:
            continue

        # 3) Probe BC (sample BC) from R2 first 8nt
        pbc_raw = seq2[:PROBE_BC_END]
        pbc = correct_probe_bc(pbc_raw)
        if pbc is None:
            continue
        n_probe_valid += 1

        # 4) Protospacer
        sgid, pos = find_protospacer(seq1, ps_to_sg, ps_1mm)
        if sgid is None:
            continue
        n_proto_valid += 1
        n_proto_pos[pos] += 1
        n_all_valid += 1

        key = (cb_corr, pbc, sgid)
        rec = counts[key]
        rec[0] += 1
        rec[1].add(umi)

        if n_total % 1_000_000 == 0:
            elapsed = time.time() - t0
            print(f"  [{time.strftime('%T')}] processed {n_total:,} reads in {elapsed:.0f}s "
                  f"(cb_valid={n_cb_valid:,}, probe_valid={n_probe_valid:,}, sgrna_valid={n_proto_valid:,})",
                  flush=True)

    fh1.close()
    fh2.close()

    print(f"[{time.strftime('%T')}] Writing output (#keys={len(counts):,})", flush=True)

    out_path = os.path.join(args.out_dir, f"{args.srr}.counts.tsv.gz")
    with gzip.open(out_path, "wt") as out:
        out.write("cell_barcode\tsample_bc\tsgrna_id\tgene\tCR\tn_reads\tn_umis\n")
        for (cb, pbc, sgid), (nr, umis) in counts.items():
            gene, cr, _ = ps_to_meta[sgid]
            out.write(f"{cb}\t{pbc}\t{sgid}\t{gene}\t{cr}\t{nr}\t{len(umis)}\n")

    stats = {
        "srr": args.srr,
        "n_total": n_total,
        "n_cb_valid": n_cb_valid,
        "n_probe_valid": n_probe_valid,
        "n_proto_valid": n_proto_valid,
        "n_all_valid": n_all_valid,
        "elapsed_s": int(time.time() - t0),
        "fraction_all_valid": n_all_valid / max(1, n_total),
        "fraction_cb_valid": n_cb_valid / max(1, n_total),
        "fraction_proto_valid": n_proto_valid / max(1, n_total),
        "proto_pos_distribution": dict(n_proto_pos),
        "n_unique_cells": len(set((k[0], k[1]) for k in counts.keys())),
        "n_unique_sgrnas_observed": len(set(k[2] for k in counts.keys())),
    }
    with open(os.path.join(args.out_dir, f"{args.srr}.stats.json"), "w") as fh:
        json.dump(stats, fh, indent=2)

    print(f"[{time.strftime('%T')}] DONE in {time.time()-t0:.0f}s")
    print(f"  total reads:     {n_total:,}")
    print(f"  CB valid:        {n_cb_valid:,}  ({100*n_cb_valid/max(1,n_total):.1f}%)")
    print(f"  probe-BC valid:  {n_probe_valid:,}  ({100*n_probe_valid/max(1,n_total):.1f}%)")
    print(f"  protospacer:     {n_proto_valid:,}  ({100*n_proto_valid/max(1,n_total):.1f}%)")
    print(f"  ALL valid:       {n_all_valid:,}  ({100*n_all_valid/max(1,n_total):.1f}%)")
    print(f"  unique (cb, sample_bc): {stats['n_unique_cells']:,}")
    print(f"  unique sgRNAs observed: {stats['n_unique_sgrnas_observed']:,}")
    print(f"  output: {out_path}")


if __name__ == "__main__":
    main()
