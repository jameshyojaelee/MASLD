#!/usr/bin/env python
"""
503_mofa_pseudobulk_prep.py — Preamble for 503_mofa_pseudobulk.R.

Builds donor × cell_type × gene pseudobulk TSVs from the atlas AnnData, one
per cell type (= one MOFA+ view). Uses HVG-subsetted counts to stay matched
with the other 5 methods in the showdown.

Output:
  results_gpu_v2/mcp/pseudobulk_mofa/view_<celltype>.tsv   (genes × donors)
  results_gpu_v2/mcp/pseudobulk_mofa/manifest.json
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
MCP = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
ATLAS_H5AD = MCP / "inputs/atlas_cnmf_global.h5ad"
HVG_FILE = MCP / "cnmf_runs/global/global.overdispersed_genes.txt"
OUT_DIR = MCP / "pseudobulk_mofa"
OUT_DIR.mkdir(parents=True, exist_ok=True)

DONOR_COL = "sample"  # per atlas metadata, 'sample' is the donor/sample id
MIN_CELLS_PER_PB = 10


def main():
    print(f"[503prep] loading {ATLAS_H5AD}")
    a = ad.read_h5ad(ATLAS_H5AD)
    print(f"[503prep] shape {a.shape}")
    hvg = set(pd.read_csv(HVG_FILE, header=None)[0].tolist())
    gene_mask = a.var_names.isin(hvg)
    print(f"[503prep] HVG genes available in atlas: {int(gene_mask.sum())}")
    # Use the counts layer (raw, for pseudobulk sum)
    counts = a.layers["counts"]
    if not sp.issparse(counts):
        counts = sp.csr_matrix(counts)
    counts = counts[:, gene_mask]
    genes = a.var_names[gene_mask]
    obs = a.obs[[DONOR_COL, "cell_type", "dataset"]].copy()
    obs["cell_type"] = obs["cell_type"].astype(str)
    obs[DONOR_COL] = obs[DONOR_COL].astype(str)

    manifest = {"views": [], "donor_col": DONOR_COL, "n_genes": int(gene_mask.sum())}
    for ct in sorted(obs["cell_type"].unique()):
        ct_idx = np.where(obs["cell_type"].values == ct)[0]
        if len(ct_idx) == 0:
            continue
        # pseudobulk per donor
        donors = obs.iloc[ct_idx][DONOR_COL].values
        uniq_donors, inv = np.unique(donors, return_inverse=True)
        # Count cells per donor; drop donors below threshold
        counts_per_donor = np.bincount(inv)
        keep = counts_per_donor >= MIN_CELLS_PER_PB
        keep_donors = uniq_donors[keep]
        if len(keep_donors) < 5:
            print(f"[503prep] view '{ct}' has <5 donors with ≥{MIN_CELLS_PER_PB} cells — skip")
            continue
        keep_mask = np.isin(donors, keep_donors)
        ct_idx_keep = ct_idx[keep_mask]
        donors_keep = donors[keep_mask]
        # pseudobulk: sum counts rows grouped by donor
        M = counts[ct_idx_keep].tocsr()
        # Build indicator matrix donor × cell
        donor_idx = pd.Categorical(donors_keep, categories=keep_donors).codes
        D = sp.csr_matrix(
            (np.ones(len(donor_idx)), (donor_idx, np.arange(len(donor_idx)))),
            shape=(len(keep_donors), M.shape[0]),
        )
        PB = (D @ M).toarray()  # donors × genes
        df = pd.DataFrame(PB.T, index=genes, columns=keep_donors)
        safe = ct.replace(" ", "_").replace("/", "_")
        out_path = OUT_DIR / f"view_{safe}.tsv"
        df.to_csv(out_path, sep="\t")
        manifest["views"].append({
            "name": ct, "path": str(out_path),
            "donor_ids": list(keep_donors), "n_genes": int(df.shape[0]),
        })
        print(f"[503prep] view '{ct}' → {df.shape[0]} genes × {df.shape[1]} donors")

    (OUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"[503prep] wrote manifest with {len(manifest['views'])} views")
    print("[503prep] DONE.")


if __name__ == "__main__":
    main()
