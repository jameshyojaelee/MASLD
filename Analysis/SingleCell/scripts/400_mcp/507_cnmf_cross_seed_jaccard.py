#!/usr/bin/env python
"""
507_cnmf_cross_seed_jaccard.py — cross-seed Jaccard noise floor.

Compares seed=42 (production) vs seed=2024 cNMF runs at k=16 by computing
top-100 Jaccard with Hungarian (linear-sum-assignment) matching between the
16 programs of run A and the 16 programs of run B.

Outputs:
  Analysis/SingleCell/results_gpu_v2/mcp/benchmarks/cnmf_cross_seed_jaccard_k16.tsv

Usage:
  python 507_cnmf_cross_seed_jaccard.py
"""
from __future__ import annotations
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
MCP = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
A = MCP / "cnmf_runs/global/global.spectra.k_16.dt_0_03.consensus.txt"
B = MCP / "cnmf_runs/global_seed2024/global_seed2024.spectra.k_16.dt_0_03.consensus.txt"
OUT = MCP / "benchmarks/cnmf_cross_seed_jaccard_k16.tsv"
TOP_N = 100


def top_genes(df: pd.DataFrame, n: int = TOP_N) -> list[set[str]]:
    """df: rows = programs, cols = genes; return top-n gene sets per program."""
    out = []
    for i in range(df.shape[0]):
        row = df.iloc[i]
        top = row.nlargest(n).index.tolist()
        out.append(set(top))
    return out


def jacc(a: set, b: set) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / max(1, len(a | b))


def main() -> None:
    if not A.exists():
        sys.exit(f"Missing seed=42 spectra: {A}")
    if not B.exists():
        sys.exit(f"Missing seed=2024 spectra: {B}")

    a = pd.read_csv(A, sep="\t", index_col=0)
    b = pd.read_csv(B, sep="\t", index_col=0)
    print(f"[507] A (seed=42)   shape: {a.shape}  programs={a.index.tolist()}")
    print(f"[507] B (seed=2024) shape: {b.shape}  programs={b.index.tolist()}")

    # Align gene namespace
    common_genes = sorted(set(a.columns) & set(b.columns))
    print(f"[507] common genes: {len(common_genes)} (A={a.shape[1]}, B={b.shape[1]})")
    a = a[common_genes]
    b = b[common_genes]

    a_top = top_genes(a)
    b_top = top_genes(b)

    # Build Jaccard matrix (16 x 16) — distance = 1 - jaccard for Hungarian min
    K = a.shape[0]
    J = np.zeros((K, K))
    for i in range(K):
        for j in range(K):
            J[i, j] = jacc(a_top[i], b_top[j])
    cost = 1.0 - J
    row_ind, col_ind = linear_sum_assignment(cost)

    matched = []
    for i, j in zip(row_ind, col_ind):
        matched.append({
            "seed42_program": a.index[i],
            "seed2024_program": b.index[j],
            "jaccard_top100": J[i, j],
            "n_intersection": len(a_top[i] & b_top[j]),
            "n_union": len(a_top[i] | b_top[j]),
        })
    df_m = pd.DataFrame(matched)
    df_m.to_csv(OUT, sep="\t", index=False, float_format="%.4f")
    print(f"[507] wrote: {OUT}")

    # Summary
    mean_j = df_m.jaccard_top100.mean()
    median_j = df_m.jaccard_top100.median()
    min_j = df_m.jaccard_top100.min()
    max_j = df_m.jaccard_top100.max()
    print(f"\n[507] cNMF cross-seed Jaccard summary (k=16, top-100 Hungarian-matched):")
    print(f"       mean   = {mean_j:.4f}")
    print(f"       median = {median_j:.4f}")
    print(f"       min    = {min_j:.4f}")
    print(f"       max    = {max_j:.4f}")
    print(f"\n[507] Interpretation: this is the noise-floor reference for the +0.040 cophenetic")
    print(f"       claim. plain-NMF-vs-cNMF Jaccard 0.89-0.92 is interpretable only relative")
    print(f"       to this floor.")
    # Append summary line for downstream
    with open(OUT, "a") as fh:
        fh.write(f"# SUMMARY\tmean={mean_j:.4f}\tmedian={median_j:.4f}\tmin={min_j:.4f}\tmax={max_j:.4f}\n")


if __name__ == "__main__":
    main()
