#!/usr/bin/env python3
"""
118a_prepare_scrna_reference.py
Prepare scRNA-seq reference for BayesPrism deconvolution.

Extracts a stratified subsample from the integrated scRNA atlas (1.23M cells)
and saves as cell-type-labeled count matrices suitable for BayesPrism.

Why subsample: BayesPrism works well with 30-50K cells; full 1.23M would be
too memory-intensive (~200GB). Stratified sampling preserves cell-type proportions.

Input:
  - Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad

Output (to results/progression/scrna_reference/):
  - scrna_counts.csv.gz      (cells × genes raw counts, subsampled)
  - scrna_cell_labels.csv     (cell_id, cell_type, condition, dataset)
  - scrna_reference_summary.csv (per-cell-type counts)

SLURM: cpu, 8 CPUs, 120G RAM, 48h
Env:   micromamba activate rapids_singlecell
"""
import os
import sys
import logging
import numpy as np
import pandas as pd
import scipy.sparse as sp

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

# ── Paths ──────────────────────────────────────────────────────────────────
BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
ATLAS = os.path.join(BASE, "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad")
OUTDIR = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression/scrna_reference")
os.makedirs(OUTDIR, exist_ok=True)

# ── Parameters ─────────────────────────────────────────────────────────────
MAX_CELLS_TOTAL = 50000       # Total cells to subsample
MAX_CELLS_PER_TYPE = 5000     # Cap per cell type to avoid hepatocyte domination
MIN_CELLS_PER_TYPE = 100      # Minimum cells per type
MIN_GENES_DETECTED = 200      # QC: minimum genes per cell
SEED = 42

# Major cell types to include (merge rare types)
CELL_TYPE_MAP = {
    "Hepatocytes": "Hepatocyte",
    "Endothelial cells": "Endothelial",
    "Fibroblasts": "Stellate",         # Stellate cells are fibroblasts in liver
    "Macrophages": "Macrophage",
    "Cholangiocytes": "Cholangiocyte",
    "T cells": "T_cell",
    "Mono+mono derived cells": "Monocyte",
    "Circulating NK/NKT": "NK_cell",
    "Resident NK": "NK_cell",          # Merge with circulating
    "Neutrophils": "Neutrophil",
    "B cells": "B_cell",
    "Plasma cells": "Plasma_cell",
    "cDC1s": "DC",
    "cDC2s": "DC",                     # Merge DC subtypes
    "pDCs": "DC",
    "Basophils": "Other_immune",
}

def main():
    log.info("=== 118a: Prepare scRNA Reference for BayesPrism ===")

    import anndata as ad
    import h5py

    # ── Load atlas metadata first (memory-efficient) ───────────────────
    log.info("Loading atlas metadata...")
    with h5py.File(ATLAS, "r") as f:
        obs = f["obs"]

        # Cell type
        ct_codes = obs["cell_type"]["codes"][:]
        ct_cats = [x.decode() if isinstance(x, bytes) else x
                   for x in obs["cell_type"]["categories"][:]]
        cell_types_raw = np.array([ct_cats[c] for c in ct_codes])

        # Condition
        cond_codes = obs["condition_harmonized"]["codes"][:]
        cond_cats = [x.decode() if isinstance(x, bytes) else x
                     for x in obs["condition_harmonized"]["categories"][:]]
        conditions = np.array([cond_cats[c] for c in cond_codes])

        # Dataset
        ds_codes = obs["dataset"]["codes"][:]
        ds_cats = [x.decode() if isinstance(x, bytes) else x
                   for x in obs["dataset"]["categories"][:]]
        datasets = np.array([ds_cats[c] for c in ds_codes])

        # Cell IDs
        cell_ids = [x.decode() if isinstance(x, bytes) else x
                    for x in obs["_index"][:]]

        n_cells = len(cell_types_raw)
        log.info(f"  Atlas: {n_cells:,} cells")

    # ── Map cell types ─────────────────────────────────────────────────
    cell_types_mapped = np.array([CELL_TYPE_MAP.get(ct, "Other") for ct in cell_types_raw])

    unique_types, type_counts = np.unique(cell_types_mapped, return_counts=True)
    log.info(f"\nMapped cell types:")
    for ct, cnt in sorted(zip(unique_types, type_counts), key=lambda x: -x[1]):
        log.info(f"  {ct}: {cnt:,}")

    # ── Stratified subsample ───────────────────────────────────────────
    log.info(f"\nStratified subsampling (max {MAX_CELLS_TOTAL:,} total, "
             f"max {MAX_CELLS_PER_TYPE:,}/type)...")

    rng = np.random.RandomState(SEED)
    selected_indices = []

    for ct in sorted(set(cell_types_mapped)):
        if ct == "Other":
            continue
        ct_mask = cell_types_mapped == ct
        ct_indices = np.where(ct_mask)[0]
        n_ct = len(ct_indices)

        if n_ct < MIN_CELLS_PER_TYPE:
            log.warning(f"  Skipping {ct}: only {n_ct} cells (< {MIN_CELLS_PER_TYPE})")
            continue

        n_sample = min(n_ct, MAX_CELLS_PER_TYPE)
        sampled = rng.choice(ct_indices, size=n_sample, replace=False)
        selected_indices.append(sampled)
        log.info(f"  {ct}: sampled {n_sample} / {n_ct}")

    selected_indices = np.concatenate(selected_indices)
    selected_indices.sort()

    # Trim to MAX_CELLS_TOTAL if needed
    if len(selected_indices) > MAX_CELLS_TOTAL:
        selected_indices = rng.choice(selected_indices, size=MAX_CELLS_TOTAL, replace=False)
        selected_indices.sort()

    log.info(f"  Total selected: {len(selected_indices):,}")

    # ── Extract raw counts for selected cells ──────────────────────────
    log.info("\nExtracting raw counts...")

    # Read the full raw.X sparse matrix
    with h5py.File(ATLAS, "r") as f:
        raw_x = f["raw"]["X"]
        data = raw_x["data"][:]
        indices = raw_x["indices"][:]
        indptr = raw_x["indptr"][:]

        # Gene names from raw.var
        raw_var = f["raw"]["var"]
        gene_names = [x.decode() if isinstance(x, bytes) else x
                      for x in raw_var["_index"][:]]

    n_genes = len(gene_names)
    log.info(f"  Sparse matrix: {n_cells:,} x {n_genes:,}")
    log.info(f"  Non-zeros: {len(data):,}")

    # Reconstruct CSR
    X_full = sp.csr_matrix((data, indices, indptr), shape=(n_cells, n_genes))

    # Subset to selected cells
    X_sub = X_full[selected_indices, :]
    log.info(f"  Subsampled matrix: {X_sub.shape}")

    # ── QC filter: minimum genes detected ──────────────────────────────
    genes_per_cell = np.diff(X_sub.indptr)  # Number of non-zero entries per row
    qc_mask = genes_per_cell >= MIN_GENES_DETECTED
    n_fail = (~qc_mask).sum()
    if n_fail > 0:
        log.info(f"  QC filtered: {n_fail} cells with < {MIN_GENES_DETECTED} genes")
        X_sub = X_sub[qc_mask, :]
        selected_indices = selected_indices[qc_mask]

    log.info(f"  Final reference: {X_sub.shape[0]:,} cells x {X_sub.shape[1]:,} genes")

    # ── Filter genes: keep genes detected in >= 10 cells ───────────────
    gene_detect = np.diff(X_sub.T.tocsr().indptr)
    gene_mask = gene_detect >= 10
    n_genes_kept = gene_mask.sum()
    log.info(f"  Genes detected in >= 10 cells: {n_genes_kept:,}")
    X_sub = X_sub[:, gene_mask]
    gene_names_after_detect = [g for g, m in zip(gene_names, gene_mask) if m]

    # ── Select top N highly variable genes (by variance) ──────────────
    # BayesPrism works well with 3-5K informative genes; 29K is too many
    MAX_GENES = 5000
    if X_sub.shape[1] > MAX_GENES:
        log.info(f"  Selecting top {MAX_GENES} highly variable genes...")
        # Compute variance per gene (on sparse matrix)
        gene_var = np.array(X_sub.power(2).mean(axis=0) - np.power(X_sub.mean(axis=0), 2)).flatten()
        top_gene_idx = np.argsort(-gene_var)[:MAX_GENES]
        top_gene_idx.sort()
        X_sub = X_sub[:, top_gene_idx]
        gene_names_filtered = [gene_names_after_detect[i] for i in top_gene_idx]
        log.info(f"  HVG selection: {X_sub.shape[1]} genes (var range: "
                 f"{gene_var[top_gene_idx].min():.2f} - {gene_var[top_gene_idx].max():.2f})")
    else:
        gene_names_filtered = gene_names_after_detect

    # ── Save cell labels ───────────────────────────────────────────────
    log.info("\nSaving results...")

    labels_df = pd.DataFrame({
        "cell_id": [cell_ids[i] for i in selected_indices],
        "cell_type": cell_types_mapped[selected_indices],
        "cell_type_raw": cell_types_raw[selected_indices],
        "condition": conditions[selected_indices],
        "dataset": datasets[selected_indices],
    })
    labels_df.to_csv(os.path.join(OUTDIR, "scrna_cell_labels.csv"), index=False)
    log.info(f"  Saved scrna_cell_labels.csv ({len(labels_df)} cells)")

    # ── Save count matrix as dense CSV (BayesPrism needs this) ─────────
    # For efficiency, save as sparse-friendly format
    log.info("  Converting to dense and saving (this may take a moment)...")
    X_dense = X_sub.toarray().astype(np.int32)

    counts_df = pd.DataFrame(
        X_dense,
        index=labels_df["cell_id"].values,
        columns=gene_names_filtered
    )
    counts_df.to_csv(os.path.join(OUTDIR, "scrna_counts.csv.gz"), compression="gzip")
    log.info(f"  Saved scrna_counts.csv.gz ({counts_df.shape})")

    # ── Summary ────────────────────────────────────────────────────────
    summary = labels_df.groupby("cell_type").size().reset_index(name="n_cells")
    summary["fraction"] = summary["n_cells"] / summary["n_cells"].sum()
    summary.to_csv(os.path.join(OUTDIR, "scrna_reference_summary.csv"), index=False)
    log.info(f"  Saved scrna_reference_summary.csv")

    log.info(f"\n=== 118a: COMPLETE ===")

if __name__ == "__main__":
    main()
