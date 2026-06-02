#!/usr/bin/env python
"""506a_compute_inverse_jaccard.py — per-cNMF-program best Jaccard from each
comparator method (plain_NMF, scHPF, LIGER).

Output rows: method, k, cnmf_program, best_jaccard_top100, best_match_factor

Used by 506b_showdown_plot.R to render Panel b heatmap and Panel d
small-multiples for the 5 paper-relevant cNMF programs (P1, P2, P11, P15, P16).

Existing `program_jaccard_*_vs_cnmf.tsv` files only record the
comparator-factor → cNMF best-match direction. This script computes the
inverse direction (cNMF program → comparator best match), which is what's
needed for a row-cNMF heatmap.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
MCP = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
BENCH = MCP / "benchmarks"
TOP_N = 100


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


def cnmf_top(k: int) -> list[tuple[str, set]]:
    p = MCP / f"cnmf_runs/global/global.spectra.k_{k}.dt_0_03.consensus.txt"
    df = pd.read_csv(p, sep="\t", index_col=0)  # rows=programs, cols=genes
    out = []
    for pidx in df.index:
        vals = df.loc[pidx].values
        order = np.argsort(vals)[::-1][:TOP_N]
        out.append((f"cnmf_P{pidx}", set(df.columns[order])))
    return out


def plain_nmf_factors(k: int) -> list[tuple[str, set]]:
    p = BENCH / f"plain_nmf_components_k{k}.npz"
    d = np.load(p, allow_pickle=True)
    comps = d["components"]  # (n_rep*k, n_genes)
    genes = d["genes"]
    out = []
    for i in range(comps.shape[0]):
        order = np.argsort(comps[i])[::-1][:TOP_N]
        rep = i // k
        f_in_rep = i % k
        out.append((f"plainNMF_rep{rep}_f{f_in_rep + 1}", set(genes[order])))
    return out


def factor_tsv_top(path: Path) -> list[tuple[str, set]]:
    df = pd.read_csv(path, sep="\t", index_col=0)  # rows=genes, cols=factors
    out = []
    for col in df.columns:
        order = df[col].values.argsort()[::-1][:TOP_N]
        out.append((col, set(df.index[order])))
    return out


def best_match_per_cnmf(
    cnmf_pp: list[tuple[str, set]],
    comparator_pp: list[tuple[str, set]],
):
    rows = []
    for cname, cset in cnmf_pp:
        best_j = 0.0
        best_f = ""
        for fname, fset in comparator_pp:
            j = jaccard(cset, fset)
            if j > best_j:
                best_j = j
                best_f = fname
        rows.append((cname, best_j, best_f))
    return rows


def main() -> int:
    out_rows = []
    for k in (10, 16):
        cnmf_pp = cnmf_top(k)

        # plain NMF
        if (BENCH / f"plain_nmf_components_k{k}.npz").exists():
            print(f"[506a] plain_NMF k={k}: loading components ...", flush=True)
            plain_pp = plain_nmf_factors(k)
            print(f"[506a] plain_NMF k={k}: {len(plain_pp)} factors; "
                  f"computing inverse Jaccard vs {len(cnmf_pp)} cNMF programs",
                  flush=True)
            for (cname, j, f) in best_match_per_cnmf(cnmf_pp, plain_pp):
                out_rows.append({
                    "method": "plain_NMF_sklearn", "k": k,
                    "cnmf_program": cname,
                    "best_jaccard_top100": j,
                    "best_match_factor": f,
                })

        # scHPF
        sf = BENCH / f"schpf_factors_k{k}.tsv"
        if sf.exists():
            print(f"[506a] scHPF k={k}: loading factors ...", flush=True)
            schpf_pp = factor_tsv_top(sf)
            for (cname, j, f) in best_match_per_cnmf(cnmf_pp, schpf_pp):
                out_rows.append({
                    "method": "scHPF", "k": k,
                    "cnmf_program": cname,
                    "best_jaccard_top100": j,
                    "best_match_factor": f,
                })

        # LIGER
        lf = BENCH / f"liger_factors_k{k}.tsv"
        if lf.exists():
            print(f"[506a] LIGER k={k}: loading factors ...", flush=True)
            liger_pp = factor_tsv_top(lf)
            for (cname, j, f) in best_match_per_cnmf(cnmf_pp, liger_pp):
                out_rows.append({
                    "method": "LIGER(pyliger)", "k": k,
                    "cnmf_program": cname,
                    "best_jaccard_top100": j,
                    "best_match_factor": f,
                })

    out_df = pd.DataFrame(out_rows)
    out_path = BENCH / "inverse_jaccard_cnmf_vs_comparators.tsv"
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"[506a] wrote {out_path} ({len(out_df)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
