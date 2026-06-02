#!/usr/bin/env python3
"""
Prepare scATAC data for Cicero co-accessibility analysis.

Strategy:
  1. Load backed h5ad (tile matrix, 500bp bins)
  2. Subset to hepatocytes (primary cell type for MASLD)
  3. Subsample to ~10K cells (Cicero recommended range)
  4. Filter to tiles accessible in >=1% of sampled cells
  5. Export as Matrix Market + metadata for R/Cicero

Also exports condition-split matrices for differential co-accessibility.

Output: Analysis/ATAC/Human_Multiome/results/cicero/
"""

import os
import sys
import logging
import numpy as np
import scipy.sparse as sp
from scipy.io import mmwrite

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
SNAP_DIR = os.path.join(BASE, "Analysis", "ATAC", "Human_Multiome", "results", "snapatac2")
OUT_DIR = os.path.join(BASE, "Analysis", "ATAC", "Human_Multiome", "results", "cicero")
os.makedirs(OUT_DIR, exist_ok=True)

N_CELLS_SAMPLE = 10000
MIN_CELL_FRAC = 0.01  # tiles accessible in >=1% of sampled cells
SEED = 42


def main():
    import anndata as ad

    logger.info("=" * 60)
    logger.info("Cicero Data Preparation")
    logger.info("=" * 60)

    # --- Load main h5ad (backed for memory efficiency) ---
    h5ad_path = os.path.join(SNAP_DIR, "snapatac2_processed.h5ad")
    logger.info(f"Loading {h5ad_path} (backed)...")
    adata = ad.read_h5ad(h5ad_path, backed="r")
    logger.info(f"  Full dataset: {adata.n_obs:,} cells x {adata.n_vars:,} features")

    # --- Get cell metadata ---
    obs = adata.obs[["cell_type", "condition", "donor_id"]].copy()
    umap = adata.obsm["X_umap"].copy()
    obs["UMAP1"] = umap[:, 0]
    obs["UMAP2"] = umap[:, 1]

    # --- Subset to hepatocytes ---
    hep_mask = obs["cell_type"] == "Hepatocyte"
    n_hep = hep_mask.sum()
    logger.info(f"  Hepatocytes: {n_hep:,} cells")

    hep_idx = np.where(hep_mask.values)[0]

    # --- Subsample cells ---
    np.random.seed(SEED)
    if n_hep > N_CELLS_SAMPLE:
        sample_idx = np.sort(np.random.choice(hep_idx, N_CELLS_SAMPLE, replace=False))
        logger.info(f"  Subsampled to {N_CELLS_SAMPLE:,} cells")
    else:
        sample_idx = hep_idx
        logger.info(f"  Using all {len(sample_idx):,} hepatocytes (< {N_CELLS_SAMPLE})")

    sample_obs = obs.iloc[sample_idx].copy()
    sample_obs.reset_index(drop=True, inplace=True)

    # --- Read tile matrix for sampled cells ---
    logger.info("Reading tile matrix for sampled cells (this may take a few minutes)...")
    # Read in chunks to manage memory
    chunk_size = 2000
    chunks = []
    for start in range(0, len(sample_idx), chunk_size):
        end = min(start + chunk_size, len(sample_idx))
        idx_slice = sample_idx[start:end]
        chunk = adata.X[idx_slice, :]
        if sp.issparse(chunk):
            chunks.append(chunk.tocsr())
        else:
            chunks.append(sp.csr_matrix(chunk))
        if (start // chunk_size) % 5 == 0:
            logger.info(f"  Read {end:,}/{len(sample_idx):,} cells...")

    X = sp.vstack(chunks, format="csr")
    logger.info(f"  Tile matrix: {X.shape[0]:,} cells x {X.shape[1]:,} tiles")

    # --- Binarize (accessibility = 0/1) ---
    X.data[:] = 1

    # --- Filter tiles by accessibility frequency ---
    min_cells = int(X.shape[0] * MIN_CELL_FRAC)
    tile_counts = np.asarray(X.sum(axis=0)).flatten()
    keep_mask = tile_counts >= min_cells
    n_keep = keep_mask.sum()
    logger.info(f"  Tiles accessible in >={MIN_CELL_FRAC*100:.0f}% of cells (>={min_cells}): {n_keep:,}")

    X_filt = X[:, keep_mask].tocsc()
    logger.info(f"  Filtered matrix: {X_filt.shape[0]:,} x {X_filt.shape[1]:,}")

    # --- Get tile coordinates ---
    var_names = list(adata.var_names)
    keep_indices = np.where(keep_mask)[0]
    tile_coords = []
    for i in keep_indices:
        name = var_names[i]  # e.g. "chr1:0-500"
        tile_coords.append(name)

    adata.file.close()

    # --- Parse tile coordinates for Cicero format ---
    # Cicero expects "chr_start_end" format
    cicero_peaks = []
    for tc in tile_coords:
        chrom, rest = tc.split(":")
        start, end = rest.split("-")
        cicero_peaks.append(f"{chrom}_{start}_{end}")

    # --- Export in Cicero triplet format ---
    # Cicero's make_atac_cds expects: peak_coord \t cell_name \t count
    logger.info("Exporting in Cicero triplet format...")

    # Binarize for export
    X_bin = X_filt.copy()
    X_bin.data[:] = 1

    cell_names = [f"cell_{i}" for i in range(len(sample_obs))]

    # Write triplet file (peak, cell, count)
    X_coo = X_bin.tocoo()
    with open(os.path.join(OUT_DIR, "cicero_input.tsv"), "w") as f:
        for i, j, v in zip(X_coo.row, X_coo.col, X_coo.data):
            f.write(f"{cicero_peaks[j]}\t{cell_names[i]}\t{int(v)}\n")
    logger.info(f"  Wrote {X_coo.nnz:,} triplets to cicero_input.tsv")

    # Peak coordinates
    with open(os.path.join(OUT_DIR, "peaks.tsv"), "w") as f:
        for p in cicero_peaks:
            f.write(p + "\n")

    # Cell barcodes
    with open(os.path.join(OUT_DIR, "barcodes.tsv"), "w") as f:
        for name in cell_names:
            f.write(name + "\n")

    # Cell metadata
    sample_obs["cell_name"] = cell_names
    sample_obs.to_csv(os.path.join(OUT_DIR, "cell_metadata.csv"), index=False)

    # --- Export per-condition triplet files for differential analysis ---
    conditions = sample_obs["condition"].unique()
    logger.info(f"Exporting per-condition files: {list(conditions)}")

    for cond in conditions:
        cond_mask = sample_obs["condition"].values == cond
        n_cond = cond_mask.sum()
        cond_safe = cond.replace(" ", "_").replace("/", "_")

        X_cond = X_bin[cond_mask, :]
        X_cond_coo = X_cond.tocoo()
        cond_cells = [cell_names[i] for i in np.where(cond_mask)[0]]

        with open(os.path.join(OUT_DIR, f"cicero_input_{cond_safe}.tsv"), "w") as f:
            for i, j, v in zip(X_cond_coo.row, X_cond_coo.col, X_cond_coo.data):
                f.write(f"{cicero_peaks[j]}\t{cond_cells[i]}\t{int(v)}\n")

        logger.info(f"  {cond}: {n_cond:,} cells")

    # --- Summary stats ---
    logger.info("\n" + "=" * 60)
    logger.info("Summary:")
    logger.info(f"  Total hepatocytes: {n_hep:,}")
    logger.info(f"  Sampled cells: {len(sample_obs):,}")
    logger.info(f"  Filtered tiles: {n_keep:,}")
    logger.info(f"  Sparsity: {1 - X_filt.nnz / (X_filt.shape[0] * X_filt.shape[1]):.4f}")
    logger.info(f"  Output: {OUT_DIR}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
