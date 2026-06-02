#!/usr/bin/env python
"""
351_pseudobulk_per_celltype_sex.py — Per-cell-type donor-level pseudobulk for sex × disease DE.

Builds donor × gene matrices for 5 cell types (hepatocyte, macrophage, fibroblast,
endothelial, cholangiocyte), summing raw counts across cells per donor. Writes
counts CSV + matched donor metadata (sample, dataset, condition, inferred_sex,
condition_binary) for downstream dream DE in Script 352.

Inputs:
  - Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad
  - Analysis/SingleCell/metadata/snrna_donor_sex.csv  (from Script 350)

Outputs (per cell type):
  - Analysis/SingleCell/results_gpu_v2/sex_celltype/{ct}_counts.csv     (genes x donors)
  - Analysis/SingleCell/results_gpu_v2/sex_celltype/{ct}_meta.csv       (donor metadata)
"""
import argparse
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

LOG_FMT = "%(asctime)s  %(levelname)s  %(message)s"
logging.basicConfig(level=logging.INFO, format=LOG_FMT, stream=sys.stderr)
log = logging.getLogger("351_pseudobulk")

CELLTYPES = {
    "Hepatocytes": "hepatocytes",
    "Macrophages": "macrophages",
    "Fibroblasts": "fibroblasts",
    "Endothelial cells": "endothelial",
    "Cholangiocytes": "cholangiocytes",
}
MIN_CELLS_PER_DONOR = 10


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--atlas", default="Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad")
    p.add_argument("--sex", default="Analysis/SingleCell/metadata/snrna_donor_sex.csv")
    p.add_argument("--outdir", default="Analysis/SingleCell/results_gpu_v2/sex_celltype")
    p.add_argument("--min-cells", type=int, default=MIN_CELLS_PER_DONOR)
    return p.parse_args()


def pseudobulk_one_celltype(adata, ct, donor_sex_df, min_cells=10):
    log.info("=== %s ===", ct)
    mask = (adata.obs["cell_type"] == ct).values
    n_cells = int(mask.sum())
    log.info("n_cells: %d", n_cells)
    if n_cells == 0:
        log.warning("no cells for %s; skipping", ct)
        return None, None

    sub = adata[mask]  # remains backed
    samples = sub.obs["sample"].astype(str).values
    uniq, inv = np.unique(samples, return_inverse=True)
    n_donors = uniq.size
    n_genes = sub.shape[1]

    counts = np.zeros((n_donors, n_genes), dtype=np.float64)
    n_per_donor = np.zeros(n_donors, dtype=np.int64)

    chunk = 100_000
    n = sub.n_obs
    for s in range(0, n, chunk):
        e = min(s + chunk, n)
        block = sub.X[s:e]
        if not sparse.issparse(block):
            block = sparse.csr_matrix(block)
        else:
            block = block.tocsr()
        inv_chunk = inv[s:e]
        # build donor x cell indicator -> aggregate
        # Approach: for each donor present in this chunk, sum rows of block belonging to it
        unique_chunk, idx = np.unique(inv_chunk, return_inverse=True)
        # Build sparse aggregation: donor_indicator (n_donors_chunk x cells_chunk)
        rows = idx
        cols = np.arange(len(inv_chunk))
        data = np.ones(len(inv_chunk), dtype=np.float64)
        agg = sparse.csr_matrix((data, (rows, cols)), shape=(len(unique_chunk), e - s))
        chunk_sums = (agg @ block).toarray()  # (n_donors_chunk, n_genes)
        for k, donor_idx in enumerate(unique_chunk):
            counts[donor_idx] += chunk_sums[k]
            n_per_donor[donor_idx] += int((idx == k).sum())
        if (s // chunk) % 3 == 0:
            log.info("  %s aggregated %d / %d cells", ct, e, n)

    counts_df = pd.DataFrame(counts.T, index=adata.var_names, columns=uniq)

    # Build meta
    obs = adata.obs[["sample", "dataset", "condition"]].astype(str).drop_duplicates(subset=["sample"])
    obs = obs.set_index("sample").loc[uniq].reset_index()
    obs["n_cells_celltype"] = n_per_donor
    obs = obs.merge(donor_sex_df[["sample", "inferred_sex", "sex_confidence"]], on="sample", how="left")

    # Condition binary mapping (matches Script 26 logic):
    #   Healthy -> Control
    #   others (MASLD, NAFLD, NASH, Cirrhotic) -> Disease
    obs["condition_binary"] = np.where(obs["condition"].isin(["Healthy"]), "Control", "Disease")

    # Filter donors below min cells
    keep = obs["n_cells_celltype"] >= min_cells
    log.info("  kept %d / %d donors with >=%d cells", int(keep.sum()), len(obs), min_cells)
    if int(keep.sum()) < 8:
        log.warning("  %s has fewer than 8 donors with sufficient cells; downstream DE may fail", ct)
    obs = obs[keep].reset_index(drop=True)
    counts_df = counts_df.loc[:, obs["sample"].tolist()]

    return counts_df, obs


def main():
    args = parse_args()
    proj_root = Path(os.environ.get("MASLD_PROJECT_ROOT",
                                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
    os.chdir(proj_root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    log.info("loading sex inference: %s", args.sex)
    sex_df = pd.read_csv(args.sex)
    log.info("sex df: %d donors, %s", len(sex_df), sex_df["inferred_sex"].value_counts().to_dict())

    log.info("opening atlas (backed): %s", args.atlas)
    adata = ad.read_h5ad(args.atlas, backed="r")

    for ct_label, ct_slug in CELLTYPES.items():
        counts_df, meta_df = pseudobulk_one_celltype(adata, ct_label, sex_df, min_cells=args.min_cells)
        if counts_df is None:
            continue
        cf = outdir / f"{ct_slug}_counts.csv.gz"
        mf = outdir / f"{ct_slug}_meta.csv"
        counts_df.to_csv(cf, compression="gzip")
        meta_df.to_csv(mf, index=False)
        log.info("wrote %s (%d genes x %d donors) and %s", cf, counts_df.shape[0], counts_df.shape[1], mf)

    log.info("done")


if __name__ == "__main__":
    main()
