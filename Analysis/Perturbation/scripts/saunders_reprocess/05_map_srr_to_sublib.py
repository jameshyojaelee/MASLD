#!/usr/bin/env python3
"""
After per-SRR demux, match each SRR to its source sublib(s) (S1_A..S8_D h5 file)
by maximizing the overlap of 24-nt cell barcodes between FASTQ-derived counts
and the h5 cell barcode set.

Each SRR contains reads from TWO of the four BC pools (BC002+BC004 or BC001+BC003,
based on empirical observation of SRR30344303/305). The mapping is:

  SRR i -> (S_w, A or C) AND (S_w, B or D)   if SRR contains BC001/BC003 pool
  SRR i -> (S_w, B or D) for BC002/BC004 pool

Each h5 file has cells with ONE specific 8-nt suffix; so we partition the SRR
cells by suffix and match each partition against h5 files separately.
"""
import argparse, glob, gzip, json, os, sys
from collections import defaultdict, Counter
import h5py


def load_h5_cells(h5_path):
    with h5py.File(h5_path, "r") as h:
        bcs = h["matrix"]["barcodes"][:]
    cells = set()
    for b in bcs:
        s = b.decode().rstrip("-1").rstrip("-")
        cells.add(s)
    return cells


def load_srr_cells(srr_counts_path):
    """Returns dict[24-nt bc] = total UMIs (across all sgRNAs)."""
    cells = defaultdict(int)
    BC_TO_RAW = {"BC001": "ACTTTAGG", "BC002": "AACGGGAA",
                 "BC003": "AGTAGGCT", "BC004": "ATGTTGAC"}
    with gzip.open(srr_counts_path, "rt") as fh:
        next(fh)
        for line in fh:
            parts = line.split("\t")
            cb16, pbc, _, _, _, _, n_umis = parts[0], parts[1], parts[2], parts[3], parts[4], int(parts[5]), int(parts[6].rstrip())
            cb24 = cb16 + BC_TO_RAW.get(pbc, pbc)
            cells[cb24] += n_umis
    return cells


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts-dir", required=True)
    ap.add_argument("--h5-glob",
                    default="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/perturbation/datasets/saunders2025/raw/h5/GSM8478327_*.h5")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    print("Loading h5 cell barcodes per sublib...")
    h5_files = sorted(glob.glob(args.h5_glob))
    h5_to_cells = {}
    for fp in h5_files:
        name = os.path.basename(fp).split("_sample")[0]
        h5_to_cells[name] = load_h5_cells(fp)
        print(f"  {name}: {len(h5_to_cells[name]):,} cells")

    print(f"\nProcessing SRR counts...")
    srr_files = sorted(glob.glob(os.path.join(args.counts_dir, "*.counts.tsv.gz")))

    rows = []
    for srr_path in srr_files:
        srr = os.path.basename(srr_path).replace(".counts.tsv.gz", "")
        srr_cells = load_srr_cells(srr_path)
        # For each h5, compute overlap
        overlaps = []
        for h5_name, h5_cells in h5_to_cells.items():
            overlap = sum(1 for c in srr_cells if c in h5_cells)
            overlaps.append((h5_name, overlap, len(h5_cells), 100*overlap/max(1,len(h5_cells))))
        overlaps.sort(key=lambda x: -x[1])
        # Print top 5
        top5 = overlaps[:5]
        print(f"\n{srr}: {len(srr_cells):,} unique cells; top 5 h5 overlaps:")
        for h5_name, ov, h5n, pct in top5:
            print(f"  {h5_name}: {ov:,}/{h5n:,} ({pct:.1f}%)")
        rows.append({
            "srr": srr,
            "n_srr_cells": len(srr_cells),
            "top1_h5": top5[0][0],
            "top1_overlap": top5[0][1],
            "top1_pct": top5[0][3],
            "top2_h5": top5[1][0],
            "top2_overlap": top5[1][1],
            "top3_h5": top5[2][0],
            "top3_overlap": top5[2][1],
        })

    with open(args.out, "w") as fh:
        json.dump(rows, fh, indent=2)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
