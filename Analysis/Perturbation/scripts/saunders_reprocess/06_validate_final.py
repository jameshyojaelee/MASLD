#!/usr/bin/env python3
"""
Validate final per-cell sgRNA assignments.

Checks:
1. % h5 cells covered (target: ≥70%)
2. NC (control) sgRNA cell fraction (expected ~11%; library is 50/456 NC)
3. sgRNA distribution (CV, max/min cell counts per sgRNA)
4. Per-CR distribution (CR1/CR2/CR3 should be ~equal)
5. Dominant fraction histogram
6. Per-sublib breakdown (S1-S8 × A-D)
"""
import argparse, gzip, glob, json, os
from collections import Counter, defaultdict
import numpy as np
import h5py


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-cell", required=True)
    ap.add_argument("--h5-glob",
                    default="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/perturbation/datasets/saunders2025/raw/h5/GSM8478327_*.h5")
    ap.add_argument("--report-out", required=True)
    args = ap.parse_args()

    # Load h5 barcodes per sublib (for per-sublib coverage)
    h5_to_cells = {}
    all_h5 = set()
    for fp in sorted(glob.glob(args.h5_glob)):
        name = os.path.basename(fp).split("_sample")[0]
        with h5py.File(fp, "r") as h:
            cells = set()
            for b in h["matrix"]["barcodes"][:]:
                s = b.decode().rstrip("-1").rstrip("-")
                cells.add(s)
            h5_to_cells[name] = cells
            all_h5 |= cells
    print(f"h5 total cells: {len(all_h5):,}; sublibs: {len(h5_to_cells)}")

    # Load assignments
    BC_TO_RAW = {"BC001": "ACTTTAGG", "BC002": "AACGGGAA",
                 "BC003": "AGTAGGCT", "BC004": "ATGTTGAC"}
    tier_count = Counter()
    sgrna_count = Counter()      # cells per sgRNA (confident only)
    gene_count = Counter()       # cells per gene
    cr_count = Counter()
    dom_fracs = []
    n_umis_list = []
    in_h5_by_tier = Counter()
    h5_by_sublib_covered = defaultdict(int)
    h5_confident_by_sublib = defaultdict(int)

    with gzip.open(args.per_cell, "rt") as fh:
        header = fh.readline().rstrip().split("\t")
        idx = {n: i for i, n in enumerate(header)}
        for line in fh:
            p = line.rstrip().split("\t")
            tier = p[idx["tier"]]
            in_h5 = p[idx["in_h5"]] == "TRUE"
            cb24 = p[idx["cell_barcode_24"]]
            sgid = p[idx["sgrna_id"]]
            gene = p[idx["target_gene"]]
            cr = p[idx["CR"]]
            n_umis = int(p[idx["n_umis"]])
            try:
                dom = float(p[idx["dominant_frac"]])
            except:
                dom = 0.0
            tier_count[tier] += 1
            in_h5_by_tier[(tier, in_h5)] += 1
            if tier == "confident":
                sgrna_count[sgid] += 1
                gene_count[gene] += 1
                cr_count[cr] += 1
                dom_fracs.append(dom)
                n_umis_list.append(n_umis)
                if in_h5:
                    # find which sublib
                    for sublib, cells in h5_to_cells.items():
                        if cb24 in cells:
                            h5_confident_by_sublib[sublib] += 1
                            break
            if in_h5:
                for sublib, cells in h5_to_cells.items():
                    if cb24 in cells:
                        h5_by_sublib_covered[sublib] += 1
                        break

    n_total = sum(tier_count.values())
    n_conf = tier_count["confident"]
    n_in_h5_conf = sum(c for (t, h), c in in_h5_by_tier.items() if t == "confident" and h)
    n_in_h5_total = sum(c for (t, h), c in in_h5_by_tier.items() if h)

    report = {
        "total_cells_with_call": n_total,
        "confident_cells": n_conf,
        "ambiguous_low_dominance": tier_count.get("ambiguous_low_dominance", 0),
        "low_total_umis": tier_count.get("low_total_umis", 0),
        "low_top_umi": tier_count.get("low_top_umi", 0),
        "h5_total_cells": len(all_h5),
        "h5_cells_with_any_call": n_in_h5_total,
        "h5_cells_confident": n_in_h5_conf,
        "pct_h5_covered_any": 100*n_in_h5_total/len(all_h5),
        "pct_h5_covered_confident": 100*n_in_h5_conf/len(all_h5),
        "unique_sgRNAs": len(sgrna_count),
        "unique_genes": len(gene_count),
        "nc_cell_count": gene_count.get("control", 0),
        "pct_nc_cells": 100*gene_count.get("control", 0)/max(1, n_conf),
        "expected_nc_pct": 100*50/456,
        "CR_distribution": dict(cr_count),
        "dominant_frac_pct50": float(np.percentile(dom_fracs, 50)) if dom_fracs else None,
        "dominant_frac_pct90": float(np.percentile(dom_fracs, 90)) if dom_fracs else None,
        "median_umis_per_cell_conf": float(np.median(n_umis_list)) if n_umis_list else None,
        "p25_umis_per_cell_conf": float(np.percentile(n_umis_list, 25)) if n_umis_list else None,
        "h5_per_sublib_coverage": {k: {"total": len(h5_to_cells[k]),
                                       "any_call": h5_by_sublib_covered[k],
                                       "confident": h5_confident_by_sublib[k],
                                       "pct_confident": 100*h5_confident_by_sublib[k]/max(1, len(h5_to_cells[k]))}
                                   for k in sorted(h5_to_cells.keys())},
        "top_15_genes": gene_count.most_common(15),
        "top_15_sgrnas": sgrna_count.most_common(15),
        "bottom_15_sgrnas": sgrna_count.most_common()[-15:],
    }

    with open(args.report_out, "w") as fh:
        json.dump(report, fh, indent=2, default=str)
    print(f"\nReport written to {args.report_out}")

    # Print summary
    print(f"\n=== VALIDATION SUMMARY ===")
    print(f"Total cells with sgRNA call: {n_total:,}")
    print(f"  Confident: {n_conf:,} ({100*n_conf/n_total:.1f}%)")
    print(f"  Ambiguous: {tier_count.get('ambiguous_low_dominance',0):,}")
    print(f"  Low UMIs:  {tier_count.get('low_total_umis',0):,}")
    print()
    print(f"h5 cell coverage: {n_in_h5_total:,}/{len(all_h5):,} ({100*n_in_h5_total/len(all_h5):.1f}%)")
    print(f"  Confident:      {n_in_h5_conf:,}/{len(all_h5):,} ({100*n_in_h5_conf/len(all_h5):.1f}%)")
    print()
    print(f"NC fraction:      {gene_count.get('control', 0)}/{n_conf} = {100*gene_count.get('control', 0)/max(1, n_conf):.1f}% (expected {100*50/456:.1f}%)")
    print(f"sgRNAs detected:  {len(sgrna_count)} / 456")
    print(f"Median UMIs/cell: {float(np.median(n_umis_list)) if n_umis_list else 0:.0f}")
    print(f"Dominant frac p50/p90: {float(np.percentile(dom_fracs, 50)) if dom_fracs else 0:.2f} / {float(np.percentile(dom_fracs, 90)) if dom_fracs else 0:.2f}")
    print()
    print(f"=== Per-sublib h5 coverage (% confident) ===")
    for k in sorted(h5_to_cells.keys()):
        cov = h5_confident_by_sublib[k]
        tot = len(h5_to_cells[k])
        print(f"  {k}: {cov:5,d}/{tot:5,d} ({100*cov/max(1,tot):.1f}%)")


if __name__ == "__main__":
    main()
