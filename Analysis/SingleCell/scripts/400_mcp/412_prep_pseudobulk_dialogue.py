#!/usr/bin/env python
"""
412_prep_pseudobulk_dialogue.py — Pseudobulk per (sample, cell_type) for DIALOGUE.

Reads the global cNMF input (410) and writes per-cell-type gene x sample raw-count TSVs
plus a combined sample metadata table. TMM + log2(CPM+1) normalization happens in R
(413_normalize_pseudobulk.R) to keep dependencies clean (edgeR is only in celltype_bio/rnaseq).

Outputs:
  results_gpu_v2/mcp/inputs/dialogue_pseudobulk/{celltype}_counts.tsv.gz
  results_gpu_v2/mcp/inputs/dialogue_pseudobulk/cell_counts_per_donor_ct.tsv
"""
from __future__ import annotations

import gzip
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
IN_ATLAS = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/atlas_cnmf_global.h5ad"
OUT_DIR = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/dialogue_pseudobulk"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MIN_CELLS_PER_DONOR_CT = 30


def pseudobulk_one_ct(adata, ct: str):
    sub = adata[adata.obs["cell_type"] == ct]
    X = sub.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    # sum counts per sample
    sample_codes, sample_uniq = pd.factorize(sub.obs["sample"], sort=True)
    n_samples = len(sample_uniq)
    n_genes = sub.n_vars
    # build row-indicator sparse matrix: samples x cells
    row = sample_codes
    col = np.arange(sub.n_obs)
    data = np.ones(sub.n_obs, dtype=np.float32)
    S = sp.csr_matrix((data, (row, col)), shape=(n_samples, sub.n_obs))
    # pseudobulk = S @ X  -> samples x genes
    P = (S @ X).toarray().astype(np.int32)
    n_cells_per_sample = np.asarray(S.sum(axis=1)).ravel().astype(int)
    # drop samples with <MIN_CELLS
    keep = n_cells_per_sample >= MIN_CELLS_PER_DONOR_CT
    P = P[keep]
    s_idx = sample_uniq[keep]
    n_cells_per_sample = n_cells_per_sample[keep]
    # genes as columns (gene_name index)
    df = pd.DataFrame(P.T, index=sub.var_names, columns=s_idx)
    df.index.name = "gene_name"
    return df, n_cells_per_sample


def main() -> None:
    print(f"[412] Loading global AnnData: {IN_ATLAS}")
    adata = ad.read_h5ad(IN_ATLAS)
    print(f"[412] Shape: {adata.shape}")

    counts_rows = []
    for ct in adata.obs["cell_type"].cat.categories if hasattr(adata.obs["cell_type"], "cat") else adata.obs["cell_type"].unique():
        if ct not in ("Hepatocytes", "Endothelial cells", "Fibroblasts", "Macrophages", "Cholangiocytes"):
            continue
        ct_safe = ct.lower().replace(" ", "_")
        print(f"[412] Pseudobulking {ct}...")
        df, n_cells = pseudobulk_one_ct(adata, ct)
        out = OUT_DIR / f"{ct_safe}_counts.tsv.gz"
        df.to_csv(out, sep="\t", compression="gzip")
        print(f"[412]   wrote {df.shape[1]} samples x {df.shape[0]} genes -> {out}")
        for s, n in zip(df.columns, n_cells):
            counts_rows.append({"sample": s, "cell_type": ct, "n_cells": int(n)})

    pd.DataFrame(counts_rows).to_csv(OUT_DIR / "cell_counts_per_donor_ct.tsv", sep="\t", index=False)
    print(f"[412] DONE. Wrote {len(counts_rows)} (donor, CT) rows.")


if __name__ == "__main__":
    main()
