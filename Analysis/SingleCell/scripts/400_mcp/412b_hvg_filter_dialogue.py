#!/usr/bin/env python
"""
412b_hvg_filter_dialogue.py — HVG pre-filter for DIALOGUE rerun.

Reads `inputs/dialogue_pseudobulk/{ct}_{counts,logcpm}.tsv.gz`, computes
per-gene CV across donors on raw counts (after CPM-normalizing per-donor to
remove library-size effects), keeps the top `N` most variable genes per cell
type, and writes `inputs/dialogue_pseudobulk_1khvg/{ct}_{counts,logcpm}.tsv.gz`
plus `inputs/dialogue_pseudobulk_1khvg/cell_counts_per_donor_ct.tsv`.

Selection per Jerby-Arnon 2022 DIALOGUE methods (~1,000 HVG per cell type).
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
IN_DIR = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/dialogue_pseudobulk"
OUT_DIR = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/dialogue_pseudobulk_1khvg"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CTS = ["hepatocytes", "endothelial_cells", "fibroblasts", "macrophages", "cholangiocytes"]
N_HVG = 1000


def select_hvg(logcpm: pd.DataFrame, n: int = N_HVG) -> list[str]:
    # logcpm: genes x donors; use CV on linear expression (mean/sd-like ranking)
    # Score genes by variance across donors in log1p-CPM space (mean-var stabilized).
    v = logcpm.var(axis=1).fillna(0.0)
    m = logcpm.mean(axis=1).fillna(0.0)
    # Keep genes with mean > 0 (expressed in at least some donors)
    expressed = m > 0.1
    score = v.where(expressed, -1)
    top = score.sort_values(ascending=False).head(n).index.tolist()
    return top


def main() -> None:
    for ct in CTS:
        counts_f = IN_DIR / f"{ct}_counts.tsv.gz"
        logcpm_f = IN_DIR / f"{ct}_logcpm.tsv.gz"
        if not (counts_f.exists() and logcpm_f.exists()):
            print(f"[412b] missing {ct}; skipping")
            continue
        counts = pd.read_csv(counts_f, sep="\t", index_col=0)
        logcpm = pd.read_csv(logcpm_f, sep="\t", index_col=0)
        print(f"[412b] {ct}: {counts.shape[0]} genes x {counts.shape[1]} donors")
        hvg = select_hvg(logcpm, N_HVG)
        counts_sub = counts.loc[counts.index.isin(hvg)]
        logcpm_sub = logcpm.loc[logcpm.index.isin(hvg)]
        counts_sub.to_csv(OUT_DIR / f"{ct}_counts.tsv.gz", sep="\t", compression="gzip")
        logcpm_sub.to_csv(OUT_DIR / f"{ct}_logcpm.tsv.gz", sep="\t", compression="gzip")
        print(f"[412b] {ct}: wrote {counts_sub.shape[0]} HVG x {counts_sub.shape[1]} donors")
    # Copy cell counts per donor/CT (unchanged)
    src = IN_DIR / "cell_counts_per_donor_ct.tsv"
    if src.exists():
        shutil.copy(src, OUT_DIR / "cell_counts_per_donor_ct.tsv")
    print("[412b] DONE.")


if __name__ == "__main__":
    main()
