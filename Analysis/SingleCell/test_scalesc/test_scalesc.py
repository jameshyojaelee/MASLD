#!/usr/bin/env python3
"""
ScaleSC GPU pipeline test using Liver Atlas Cell Ranger outputs.

Steps:
  1. Convert Cell Ranger filtered_feature_bc_matrix.h5 → per-sample h5ad files
  2. Run ScaleSC full pipeline (QC → HVG → normalize → PCA → Harmony → neighbors → leiden → UMAP)
  3. Run cluster merging + marker identification
  4. Compare timing vs rapids_singlecell baseline

Usage:
  python test_scalesc.py --output-dir /path/to/output [--n-samples 5] [--max-cell-batch 100000]
"""

import argparse
import logging
import os
import sys
import time
from pathlib import Path

import numpy as np
import scanpy as sc
import anndata as ad

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def find_cellranger_h5(base_dir: str, n_samples: int | None = None) -> list[tuple[str, str]]:
    """Find Cell Ranger filtered h5 files, return list of (sample_id, h5_path)."""
    base = Path(base_dir)
    samples = []
    for h5 in sorted(base.glob("*/outs/filtered_feature_bc_matrix.h5")):
        sample_id = h5.parent.parent.name
        samples.append((sample_id, str(h5)))
    if n_samples:
        samples = samples[:n_samples]
    logger.info(f"Found {len(samples)} Cell Ranger samples")
    return samples


def prepare_h5ad_dir(samples: list[tuple[str, str]], output_dir: str) -> str:
    """Convert Cell Ranger h5 files to per-sample h5ad files in a directory for ScaleSC."""
    h5ad_dir = os.path.join(output_dir, "scalesc_input")
    os.makedirs(h5ad_dir, exist_ok=True)

    total_cells = 0
    for sample_id, h5_path in samples:
        out_path = os.path.join(h5ad_dir, f"{sample_id}.h5ad")
        if os.path.exists(out_path):
            adata = sc.read_h5ad(out_path)
            logger.info(f"  {sample_id}: {adata.n_obs} cells (cached)")
            total_cells += adata.n_obs
            continue

        t0 = time.time()
        adata = sc.read_10x_h5(h5_path)
        adata.var_names_make_unique()
        adata.obs["sample"] = sample_id
        # Ensure sparse CSR format (ScaleSC expects this)
        import scipy.sparse
        if not scipy.sparse.issparse(adata.X):
            adata.X = scipy.sparse.csr_matrix(adata.X)
        elif not scipy.sparse.isspmatrix_csr(adata.X):
            adata.X = adata.X.tocsr()

        adata.write_h5ad(out_path)
        logger.info(f"  {sample_id}: {adata.n_obs} cells ({time.time()-t0:.1f}s)")
        total_cells += adata.n_obs

    logger.info(f"Total: {total_cells} cells across {len(samples)} samples → {h5ad_dir}")
    return h5ad_dir


def run_scalesc_pipeline(h5ad_dir: str, output_dir: str, max_cell_batch: int = 100000):
    """Run the full ScaleSC GPU pipeline."""
    import scalesc as ssc

    results_dir = os.path.join(output_dir, "scalesc_results")
    os.makedirs(results_dir, exist_ok=True)

    logger.info("=" * 60)
    logger.info("ScaleSC Pipeline Start")
    logger.info("=" * 60)
    t_total = time.time()

    # Initialize ScaleSC
    t0 = time.time()
    pipeline = ssc.ScaleSC(
        data_dir=h5ad_dir,
        max_cell_batch=int(max_cell_batch),
        preload_on_cpu=True,
        preload_on_gpu=True,
        save_raw_counts=False,
        save_norm_counts=False,
        save_after_each_step=False,
        output_dir=results_dir,
    )
    logger.info(f"Init + data load: {time.time()-t0:.1f}s")

    # QC
    t0 = time.time()
    pipeline.calculate_qc_metrics()
    logger.info(f"QC metrics: {time.time()-t0:.1f}s | Shape: {pipeline.adata.shape}")

    # Filter genes and cells
    t0 = time.time()
    pipeline.filter_genes(min_count=3)
    pipeline.filter_cells(min_count=200, max_count=6000)
    logger.info(f"Filtering: {time.time()-t0:.1f}s | Shape after: {pipeline.adata.shape}")

    # HVG
    t0 = time.time()
    pipeline.highly_variable_genes(n_top_genes=4000)
    n_hvg = pipeline.adata.var["highly_variable"].sum()
    logger.info(f"HVG selection: {time.time()-t0:.1f}s | {n_hvg} HVGs")

    # Normalize + log1p
    t0 = time.time()
    pipeline.normalize_log1p()
    logger.info(f"Normalize + log1p: {time.time()-t0:.1f}s")

    # PCA
    t0 = time.time()
    pipeline.pca(n_components=50)
    logger.info(f"PCA: {time.time()-t0:.1f}s")

    # Harmony batch correction
    t0 = time.time()
    pipeline.harmony(sample_col_name="sample", max_iter_harmony=20)
    logger.info(f"Harmony: {time.time()-t0:.1f}s")

    # Neighbors
    t0 = time.time()
    pipeline.neighbors(n_neighbors=20, n_pcs=50, use_rep="X_pca_harmony")
    logger.info(f"Neighbors: {time.time()-t0:.1f}s")

    # Leiden clustering
    t0 = time.time()
    pipeline.leiden(resolution=0.5)
    n_clusters = pipeline.adata.obs["leiden"].nunique()
    logger.info(f"Leiden: {time.time()-t0:.1f}s | {n_clusters} clusters")

    # UMAP
    t0 = time.time()
    pipeline.umap()
    logger.info(f"UMAP: {time.time()-t0:.1f}s")

    # Save results
    t0 = time.time()
    pipeline.save(data_name="liver_atlas_scalesc")
    logger.info(f"Save: {time.time()-t0:.1f}s")

    total_time = time.time() - t_total
    logger.info("=" * 60)
    logger.info(f"ScaleSC Pipeline Complete: {total_time:.1f}s total")
    logger.info(f"  Cells: {pipeline.adata.n_obs}, Genes: {pipeline.adata.n_vars}")
    logger.info(f"  Clusters: {n_clusters}")
    logger.info("=" * 60)

    # --- Cluster merging + marker identification ---
    logger.info("\n--- Cluster Merging + Marker Identification ---")
    t0 = time.time()
    pipeline.to_CPU()
    adata_full = pipeline.adata_X
    adata_full.obs["leiden"] = pipeline.adata.obs["leiden"]
    adata_full.obsm["X_umap"] = pipeline.adata.obsm["X_umap"]
    adata_full.var = pipeline.adata.var
    adata_full = adata_full[:, adata_full.var["highly_variable"]].copy()

    merged_col = ssc.clusters_merge(adata_full, "leiden")
    n_merged = adata_full.obs[merged_col].nunique()
    logger.info(f"Cluster merging: {time.time()-t0:.1f}s | {n_clusters} → {n_merged} clusters")

    # Save the full annotated object
    adata_out = pipeline.adata.copy()
    adata_out.obs[merged_col] = adata_full.obs[merged_col]
    out_path = os.path.join(results_dir, "liver_atlas_scalesc_annotated.h5ad")
    adata_out.write_h5ad(out_path)
    logger.info(f"Saved annotated atlas to {out_path}")

    return total_time


def run_rapids_baseline(h5ad_dir: str, output_dir: str):
    """Run rapids_singlecell baseline for comparison."""
    import rapids_singlecell as rsc
    from rapids_singlecell.get import anndata_to_GPU, anndata_to_CPU

    logger.info("=" * 60)
    logger.info("rapids_singlecell Baseline Start")
    logger.info("=" * 60)
    t_total = time.time()

    # Load all h5ad files
    t0 = time.time()
    adatas = []
    for f in sorted(Path(h5ad_dir).glob("*.h5ad")):
        a = sc.read_h5ad(str(f))
        adatas.append(a)
    adata = ad.concat(adatas, merge="same")
    adata.var_names_make_unique()
    logger.info(f"Load + concat: {time.time()-t0:.1f}s | {adata.shape}")

    # QC + filter
    t0 = time.time()
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_genes=200)
    adata = adata[adata.obs["n_genes"] < 6000].copy()
    logger.info(f"Filter: {time.time()-t0:.1f}s | {adata.shape}")

    # HVG on CPU (seurat_v3 requires counts)
    t0 = time.time()
    sc.pp.highly_variable_genes(adata, n_top_genes=4000, flavor="seurat_v3")
    logger.info(f"HVG: {time.time()-t0:.1f}s")

    # Normalize
    t0 = time.time()
    anndata_to_GPU(adata)
    rsc.pp.normalize_total(adata, target_sum=1e4)
    rsc.pp.log1p(adata)
    logger.info(f"Normalize (GPU): {time.time()-t0:.1f}s")

    # PCA + neighbors + UMAP + leiden (all GPU)
    t0 = time.time()
    rsc.pp.scale(adata)
    rsc.pp.pca(adata, n_comps=50)
    logger.info(f"Scale + PCA (GPU): {time.time()-t0:.1f}s")

    t0 = time.time()
    rsc.pp.neighbors(adata, n_neighbors=20, n_pcs=50)
    logger.info(f"Neighbors (GPU): {time.time()-t0:.1f}s")

    t0 = time.time()
    rsc.tl.leiden(adata, resolution=0.5)
    logger.info(f"Leiden (GPU): {time.time()-t0:.1f}s")

    t0 = time.time()
    rsc.tl.umap(adata)
    logger.info(f"UMAP (GPU): {time.time()-t0:.1f}s")

    anndata_to_CPU(adata)
    total_time = time.time() - t_total
    n_clusters = adata.obs["leiden"].nunique()

    logger.info("=" * 60)
    logger.info(f"rapids_singlecell Baseline Complete: {total_time:.1f}s total")
    logger.info(f"  Cells: {adata.n_obs}, Clusters: {n_clusters}")
    logger.info("=" * 60)

    return total_time


def main():
    parser = argparse.ArgumentParser(description="ScaleSC GPU pipeline test")
    parser.add_argument("--output-dir", required=True, help="Output directory for results")
    parser.add_argument("--cellranger-dir",
                        default="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Liver_Atlas/cellranger/human",
                        help="Cell Ranger output directory")
    parser.add_argument("--n-samples", type=int, default=None, help="Limit to N samples (default: all)")
    parser.add_argument("--max-cell-batch", type=int, default=100000, help="Max cells per ScaleSC batch")
    parser.add_argument("--skip-baseline", action="store_true", help="Skip rapids_singlecell baseline")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Step 1: Find and convert Cell Ranger outputs
    logger.info("Step 1: Preparing h5ad files from Cell Ranger outputs")
    samples = find_cellranger_h5(args.cellranger_dir, args.n_samples)
    h5ad_dir = prepare_h5ad_dir(samples, args.output_dir)

    # Step 2: Run ScaleSC pipeline
    logger.info("\nStep 2: Running ScaleSC pipeline")
    scalesc_time = run_scalesc_pipeline(h5ad_dir, args.output_dir, args.max_cell_batch)

    # Step 3: Run rapids_singlecell baseline
    if not args.skip_baseline:
        logger.info("\nStep 3: Running rapids_singlecell baseline")
        rapids_time = run_rapids_baseline(h5ad_dir, args.output_dir)

        logger.info("\n" + "=" * 60)
        logger.info("COMPARISON")
        logger.info(f"  ScaleSC:            {scalesc_time:.1f}s")
        logger.info(f"  rapids_singlecell:  {rapids_time:.1f}s")
        logger.info(f"  Ratio:              {rapids_time/scalesc_time:.2f}x")
        logger.info("=" * 60)
    else:
        logger.info(f"\nScaleSC completed in {scalesc_time:.1f}s")


if __name__ == "__main__":
    main()
