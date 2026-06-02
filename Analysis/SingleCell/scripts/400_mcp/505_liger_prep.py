#!/usr/bin/env python
"""
505_liger_prep.py — Stage per-dataset HVG count matrices for LIGER.

Splits atlas_cnmf_global.h5ad by `dataset` (7 datasets), subsets to HVG
(2,499 overdispersed genes), writes one .mtx.gz + _genes.txt + _barcodes.txt
per dataset plus a manifest.json.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
MCP = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
ATLAS = MCP / "inputs/atlas_cnmf_global.h5ad"
HVG = MCP / "cnmf_runs/global/global.overdispersed_genes.txt"
OUT = MCP / "liger_staging"
OUT.mkdir(parents=True, exist_ok=True)

SUB_PER_DATASET = int(os.environ.get("LIGER_N_CELLS_PER_DATASET", 0))  # 0 = full


def main():
    print(f"[505prep] loading {ATLAS}")
    a = ad.read_h5ad(ATLAS)
    hvg_set = set(pd.read_csv(HVG, header=None)[0].tolist())
    mask = a.var_names.isin(hvg_set)
    X = a.layers["counts"][:, mask]
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    genes = a.var_names[mask]
    obs = a.obs.copy()
    obs["dataset"] = obs["dataset"].astype(str)
    datasets = sorted(obs["dataset"].unique())
    print(f"[505prep] datasets: {datasets}")

    rng = np.random.default_rng(42)
    manifest = {"datasets": [], "n_genes": int(mask.sum())}
    for d in datasets:
        idx = np.where(obs["dataset"].values == d)[0]
        if SUB_PER_DATASET and len(idx) > SUB_PER_DATASET:
            idx = rng.choice(idx, SUB_PER_DATASET, replace=False)
            idx.sort()
        sub = X[idx].T.tocoo()  # genes × cells
        barcodes = a.obs_names[idx].tolist()
        tag = d.replace("/", "_").replace(" ", "_")
        mtx_path = OUT / f"{tag}.mtx"
        mmwrite(str(mtx_path), sub)
        # gzip to save space — but rliger can read plain .mtx too; keep plain
        genes_path = OUT / f"{tag}_genes.txt"
        genes_path.write_text("\n".join(genes.tolist()) + "\n")
        barcodes_path = OUT / f"{tag}_barcodes.txt"
        barcodes_path.write_text("\n".join(str(b) for b in barcodes) + "\n")
        manifest["datasets"].append({
            "name": d,
            "mtx_path": str(mtx_path),
            "genes_path": str(genes_path),
            "barcodes_path": str(barcodes_path),
            "n_cells": int(len(idx)),
        })
        print(f"[505prep]   {d}: {sub.shape[0]}×{sub.shape[1]} ({len(idx)} cells)")

    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"[505prep] wrote manifest with {len(datasets)} datasets")
    print("[505prep] DONE.")


if __name__ == "__main__":
    main()
