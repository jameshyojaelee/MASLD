#!/usr/bin/env python3
"""
300: Pseudotime Preprocessing — Subset & Re-embed 5 Core Cell Types.

Subsets the integrated scRNA-seq atlas (1.23M cells) into 5 cell-type-specific
objects, each re-normalized, re-embedded (PCA/UMAP), and re-clustered at higher
resolution to identify disease sub-states for trajectory analysis.

Cell types: Hepatocytes, Macrophages, Fibroblasts, Endothelial cells, Cholangiocytes

Inputs:
    - Integrated atlas: Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad

Outputs (to results_gpu_v2/pseudotime/):
    - {celltype}_subset.h5ad (5 files)
    - preprocessing_summary.csv

Usage:
    sbatch run_pseudotime_pipeline.sh  # or run this script directly on GPU node
"""

import os
import sys
import warnings
import logging
import gc
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
SC_DIR = os.path.join(BASE, "Analysis/SingleCell")
RESULTS = os.path.join(SC_DIR, "results_gpu_v2")
H5AD = os.path.join(RESULTS, "integrated_atlas.h5ad")
OUT_DIR = os.path.join(RESULTS, "pseudotime")
os.makedirs(OUT_DIR, exist_ok=True)

# Target cell types
CELL_TYPES = [
    "Hepatocytes",
    "Macrophages",
    "Fibroblasts",
    "Endothelial cells",
    "Cholangiocytes",
]

# Subsampling threshold for large cell types
MAX_CELLS = 100_000

# ---------------------------------------------------------------------------
# GPU init
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.join(SC_DIR, "scripts"))
try:
    from gpu_utils import init_gpu, get_processor
    USE_GPU = init_gpu()
except Exception:
    USE_GPU = False
    log.info("No GPU — running CPU-only")

# ---------------------------------------------------------------------------
# Imports (after GPU init)
# ---------------------------------------------------------------------------
import scanpy as sc
import anndata as ad
from scipy import sparse

if USE_GPU:
    pp, tl = get_processor(True)
    from gpu_utils import to_gpu, from_gpu
else:
    pp, tl = sc.pp, sc.tl

# ---------------------------------------------------------------------------
# Fix h5ad categorical encoding (missing 'ordered' attribute)
# ---------------------------------------------------------------------------
def read_h5ad_fixed(h5ad_path):
    """Read h5ad with categorical 'ordered' attribute fix (in-memory, no file modification)."""
    import h5py as _h5py
    # First, check if fix is needed by trying a direct read
    try:
        return sc.read_h5ad(h5ad_path)
    except KeyError:
        pass

    # Fix: copy file to temp, patch, then read
    import shutil
    import tempfile
    tmp_path = os.path.join(OUT_DIR, ".atlas_patched.h5ad")
    if not os.path.exists(tmp_path):
        log.info("Patching h5ad categoricals (one-time copy to %s)...", tmp_path)
        shutil.copy2(h5ad_path, tmp_path)
        with _h5py.File(tmp_path, "a") as f:
            for col in f["obs"].keys():
                if col == "_index":
                    continue
                attrs = f["obs"][col].attrs
                enc = attrs.get("encoding-type", b"")
                if isinstance(enc, bytes):
                    enc = enc.decode()
                if enc == "categorical" and "ordered" not in attrs:
                    attrs["ordered"] = False
                    log.info("  Fixed '%s'", col)
    return sc.read_h5ad(tmp_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    log.info("=" * 60)
    log.info("300: Pseudotime Preprocessing")
    log.info("=" * 60)

    # Load atlas (with categorical fix — does NOT modify original file)
    log.info("Loading atlas from %s ...", H5AD)
    adata = read_h5ad_fixed(H5AD)
    log.info("Atlas shape: %s", adata.shape)
    log.info("Cell types: %s", adata.obs["cell_type"].value_counts().to_dict())
    log.info("Conditions: %s", adata.obs["condition"].value_counts().to_dict())

    # Extract all 5 subsets in one pass, then release the full atlas
    raw_subsets = {}
    for ct in CELL_TYPES:
        mask = adata.obs["cell_type"] == ct
        n_total = mask.sum()
        log.info("%s: %d cells", ct, n_total)
        if n_total > 0:
            raw_subsets[ct] = (adata[mask].copy(), n_total)

    # FREE the full atlas from memory
    del adata
    gc.collect()
    log.info("Released atlas from memory — processing %d subsets", len(raw_subsets))

    summary_rows = []

    for ct in CELL_TYPES:
        log.info("\n" + "=" * 50)
        log.info("Processing: %s", ct)
        log.info("=" * 50)

        if ct not in raw_subsets:
            log.warning("No %s cells — skipping", ct)
            continue

        sub, n_total = raw_subsets.pop(ct)  # pop to free memory after processing
        log.info("Total %s cells: %d", ct, n_total)

        # Filter low-confidence cells if confidence column exists
        if "cell_type_conf" in sub.obs.columns:
            high_conf = sub.obs["cell_type_conf"].astype(float) > 0.5
            n_low = (~high_conf).sum()
            if n_low > 0:
                log.info("Removing %d low-confidence cells (%.1f%%)",
                         n_low, 100 * n_low / len(sub))
                sub = sub[high_conf].copy()

        # Subsample if too large
        n_cells = len(sub)
        subsampled = False
        if n_cells > MAX_CELLS:
            log.info("Subsampling from %d to %d cells", n_cells, MAX_CELLS)
            np.random.seed(42)
            idx = np.random.choice(n_cells, MAX_CELLS, replace=False)
            idx.sort()
            sub = sub[idx].copy()
            subsampled = True
            n_cells = len(sub)

        log.info("Working with %d cells", n_cells)

        # NOTE: Atlas X is already log-normalized. Do NOT re-normalize.
        # Use flavor="seurat" which works on log-normalized data (not "seurat_v3"
        # which expects raw counts in a layer).
        log.info("Data is already log-normalized (atlas X range: [%.2f, %.2f])",
                 sub.X.min() if not sparse.issparse(sub.X) else sub.X.data.min(),
                 sub.X.max() if not sparse.issparse(sub.X) else sub.X.data.max())

        if USE_GPU:
            to_gpu(sub)

        # HVG selection (within cell type) — seurat flavor on log-normalized data
        log.info("Selecting HVGs...")
        n_hvg = min(3000, sub.shape[1] - 1)
        sc.pp.highly_variable_genes(sub, n_top_genes=n_hvg, flavor="seurat")
        log.info("HVGs selected: %d", sub.var["highly_variable"].sum())

        # Scale
        pp.scale(sub, max_value=10)

        # PCA
        log.info("Computing PCA...")
        tl.pca(sub, n_comps=50)

        # Neighbors + UMAP
        log.info("Computing neighbors + UMAP...")
        pp.neighbors(sub, n_neighbors=30, n_pcs=50)
        tl.umap(sub)

        # Re-cluster at higher resolution for sub-state identification
        log.info("Clustering (resolution=1.0)...")
        tl.leiden(sub, resolution=1.0, key_added="leiden_substate")

        if USE_GPU:
            from_gpu(sub)

        n_clusters = sub.obs["leiden_substate"].nunique()
        log.info("Found %d sub-clusters", n_clusters)

        # Log condition distribution
        cond_dist = sub.obs["condition"].value_counts().to_dict()
        log.info("Condition distribution: %s", cond_dist)

        # Save h5ad
        ct_safe = ct.replace(" ", "_").replace("/", "_")
        out_path = os.path.join(OUT_DIR, f"{ct_safe}_subset.h5ad")
        log.info("Saving to %s ...", out_path)
        sub.write_h5ad(out_path)

        # Export metadata + UMAP for R consumption (avoids reticulate dependency)
        meta_df = sub.obs.copy()
        if "X_umap" in sub.obsm:
            meta_df["UMAP_1"] = sub.obsm["X_umap"][:, 0]
            meta_df["UMAP_2"] = sub.obsm["X_umap"][:, 1]
        meta_path = os.path.join(OUT_DIR, f"{ct_safe}_metadata.csv")
        meta_df.to_csv(meta_path)
        log.info("Exported metadata+UMAP CSV for R: %s", meta_path)

        # NOTE: h5ad files are used directly by R via reticulate (no MTX export
        # needed — text-based MatrixMarket is impractical for 37K x 100K matrices)

        summary_rows.append({
            "cell_type": ct,
            "n_total": n_total,
            "n_after_qc": n_cells,
            "subsampled": subsampled,
            "n_clusters": n_clusters,
            **{f"n_{k}": v for k, v in cond_dist.items()},
        })

        del sub
        gc.collect()

    # Save summary
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(os.path.join(OUT_DIR, "preprocessing_summary.csv"), index=False)
    log.info("\nPreprocessing summary:\n%s", summary.to_string())

    log.info("\n=== 300: Preprocessing COMPLETE ===")


if __name__ == "__main__":
    main()
