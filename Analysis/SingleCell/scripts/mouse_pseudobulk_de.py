#!/usr/bin/env python3
"""
Mouse scRNA-seq Pseudobulk Generation for Species-Effect-Size Comparison.

This script:
1. Loads per-sample mouse h5ad files (raw counts)
2. Annotates cells with liver cell types using marker gene scoring
3. Aggregates raw counts per (sample, cell_type) → pseudobulk CSVs

Context:
  The integrated mouse h5ad (scalesc_mouse_annotated.h5ad) has leiden clusters
  but no saved counts (X=None). Raw counts are in the per-sample input h5ads.
  Marker-based annotation is applied fresh to each sample.

Design notes:
  - Only 4 GSE189600 samples have clear MASLD-like (ALIOS) vs Healthy
    (Normal_Chow) condition labels.
  - Liver_Atlas mouse samples (N=110) are labeled "Mixed" — excluded from
    binary DE but included in pseudobulk aggregation for reference.
  - Pseudobulk CSVs are saved for all samples; the R DE script filters to
    GSE189600 (ALIOS vs Normal_Chow) for the actual DE analysis.

Output directory: results_gpu_v2/mouse_sc/
  pseudobulk/{celltype}_pseudobulk.csv  — genes x samples (raw counts)

Usage (via SLURM):
  python mouse_pseudobulk_de.py
"""

import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SC_DIR = BASE / "Analysis/SingleCell"
INPUT_DIR = SC_DIR / "integration/input_h5ad/mouse"
OUT_DIR = SC_DIR / "results_gpu_v2/mouse_sc/pseudobulk"
MANIFEST = SC_DIR / "integration/input_h5ad/sample_manifest.csv"

OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Mouse liver cell type marker gene sets
# ---------------------------------------------------------------------------
MOUSE_MARKERS = {
    "Hepatocytes": [
        "Alb", "Apoe", "Cyp2e1", "Cyp3a11", "Tat", "Arg1", "Ass1",
        "Serpina1a", "Ttr", "Apob", "Fabp1", "Apoa1", "Cyp7a1", "Hmgcr",
    ],
    "Cholangiocytes": [
        "Krt7", "Krt19", "Epcam", "Spp1", "Anxa4", "Sox9", "Cftr",
    ],
    "LSECs": [
        "Pecam1", "Cd34", "Vwf", "Lyve1", "Stab2", "Oit3", "Fcgr2b",
    ],
    "Kupffer_cells": [
        "Csf1r", "Adgre1", "Clec4f", "Timd4", "Vsig4", "Cd68", "Cx3cr1",
    ],
    "Stellate_cells": [
        "Acta2", "Col1a1", "Dcn", "Pdgfrb", "Lrat", "Des", "Reln",
    ],
    "T_cells": [
        "Cd3d", "Cd3e", "Trdc", "Cd8a", "Cd4", "Nkg7",
    ],
    "NK_cells": [
        "Nkg7", "Gzma", "Gzmb", "Ncr1", "Klrb1c",
    ],
    "B_cells": [
        "Cd79a", "Ms4a1", "Cd19", "Ighm",
    ],
    "Monocytes": [
        "Ly6c2", "Ccr2", "Csf3r", "S100a8", "S100a9",
    ],
    "pDCs": [
        "Siglech", "Bst2", "Irf7",
    ],
}


def score_cell_types(counts_csr: sp.csr_matrix, gene_names: list) -> pd.DataFrame:
    """
    Score each cell for each marker gene set using mean normalized expression.
    Returns a DataFrame (cells x cell_types) with scores.

    Method:
      1. Normalize counts to 10K (median), log1p
      2. For each cell type, compute mean expression of available markers
    """
    # Normalize: 10K normalization + log1p (fast, in-place on dense slice)
    counts_csr = counts_csr.astype(np.float32)
    row_sums = np.asarray(counts_csr.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1.0  # avoid division by zero
    # Scale factor per cell
    scale = 10000.0 / row_sums

    # Apply normalization row-wise (efficiently for sparse)
    from scipy.sparse import diags
    norm_mat = diags(scale) @ counts_csr
    # log1p
    norm_mat.data = np.log1p(norm_mat.data)

    gene_index = {g: i for i, g in enumerate(gene_names)}
    n_cells = counts_csr.shape[0]

    scores = {}
    for ct, markers in MOUSE_MARKERS.items():
        avail = [g for g in markers if g in gene_index]
        if not avail:
            scores[ct] = np.zeros(n_cells, dtype=np.float32)
            continue
        idx = [gene_index[g] for g in avail]
        # Mean across available markers (dense extraction for subset)
        sub = norm_mat[:, idx].toarray()
        scores[ct] = sub.mean(axis=1).astype(np.float32)

    return pd.DataFrame(scores)


def assign_cell_types(score_df: pd.DataFrame) -> np.ndarray:
    """Assign each cell to the highest-scoring cell type."""
    return score_df.idxmax(axis=1).values


def process_sample(h5ad_path: Path, sample_id: str) -> dict | None:
    """
    Load one sample h5ad, annotate cells, return dict of
    {cell_type: 1D count array (summed over assigned cells)} and gene_names.
    """
    import anndata as ad

    try:
        adata = ad.read_h5ad(str(h5ad_path))
    except Exception as e:
        log.warning(f"  Could not read {h5ad_path}: {e}")
        return None

    if adata.X is None or adata.shape[0] == 0:
        log.warning(f"  Empty or missing counts in {h5ad_path}")
        return None

    n_cells, n_genes = adata.shape
    gene_names = list(adata.var_names)

    # Get sparse counts
    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    elif not isinstance(X, sp.csr_matrix):
        X = X.tocsr()

    log.info(f"  {sample_id}: {n_cells:,} cells x {n_genes:,} genes")

    # Score and assign cell types
    score_df = score_cell_types(X, gene_names)
    assignments = assign_cell_types(score_df)

    # Aggregate raw counts per cell type
    result = {"gene_names": gene_names}
    for ct in MOUSE_MARKERS:
        mask = assignments == ct
        n_ct = mask.sum()
        if n_ct > 0:
            agg = np.asarray(X[mask].sum(axis=0)).ravel()
        else:
            agg = np.zeros(n_genes, dtype=np.float64)
        result[ct] = agg

    ct_counts = {ct: (assignments == ct).sum() for ct in MOUSE_MARKERS}
    log.info(f"  Cell type distribution: "
             f"{', '.join(f'{k}={v}' for k, v in sorted(ct_counts.items(), key=lambda x: -x[1])[:5])}")

    return result


def main():
    log.info("=== Mouse scRNA-seq Pseudobulk Generation ===")

    # Load manifest — filter to mouse samples
    manifest = pd.read_csv(MANIFEST)
    manifest_mouse = manifest[manifest["species"] == "mouse"].copy()
    log.info(f"Mouse samples in manifest: {len(manifest_mouse)}")
    log.info(f"Conditions: {manifest_mouse['condition'].value_counts().to_dict()}")

    # Match manifest to available h5ad files
    available = {p.stem: p for p in INPUT_DIR.glob("*.h5ad")}
    log.info(f"Available h5ad files: {len(available)}")

    manifest_mouse = manifest_mouse[manifest_mouse["srr"].isin(available)].copy()
    log.info(f"Matched: {len(manifest_mouse)} samples")

    # Process each sample
    all_results = {}  # sample_id -> {cell_type -> count_array}
    gene_names_ref = None

    for _, row in manifest_mouse.iterrows():
        srr = row["srr"]
        condition = row["condition"]
        dataset = row["dataset"]
        h5ad_path = available[srr]

        log.info(f"Processing {srr} (dataset={dataset}, condition={condition})")
        result = process_sample(h5ad_path, srr)

        if result is None:
            continue

        if gene_names_ref is None:
            gene_names_ref = result["gene_names"]

        # Check gene consistency
        if result["gene_names"] != gene_names_ref:
            log.warning(f"  Gene names differ for {srr} — skipping")
            continue

        all_results[srr] = {
            "condition": condition,
            "dataset": dataset,
            **{ct: result[ct] for ct in MOUSE_MARKERS}
        }

    if not all_results:
        log.error("No samples processed successfully.")
        sys.exit(1)

    log.info(f"\nSuccessfully processed: {len(all_results)} samples")

    # Save per-cell-type pseudobulk matrices (genes x samples)
    sample_ids = list(all_results.keys())
    conditions = [all_results[s]["condition"] for s in sample_ids]
    datasets = [all_results[s]["dataset"] for s in sample_ids]

    for ct in MOUSE_MARKERS:
        # Build count matrix: genes x samples
        count_arrays = [all_results[s][ct] for s in sample_ids]
        count_mat = np.column_stack(count_arrays)  # genes x samples

        # Create DataFrame
        df = pd.DataFrame(count_mat, index=gene_names_ref, columns=sample_ids)
        df.index.name = "gene"

        out_path = OUT_DIR / f"{ct}_pseudobulk.csv"
        df.to_csv(out_path)
        log.info(f"  Saved: {out_path} ({df.shape[0]:,} genes x {df.shape[1]} samples)")

    # Save sample metadata
    meta_df = pd.DataFrame({
        "sample": sample_ids,
        "condition": conditions,
        "dataset": datasets,
    })
    meta_path = OUT_DIR.parent / "sample_metadata.csv"
    meta_df.to_csv(meta_path, index=False)
    log.info(f"  Saved: {meta_path}")

    log.info("=== Pseudobulk generation complete ===")
    log.info(f"Output directory: {OUT_DIR}")


if __name__ == "__main__":
    main()
