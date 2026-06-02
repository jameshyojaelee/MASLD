#!/usr/bin/env python3
"""
308: Hepatocyte Subclustering via scVI

Train scVI on hepatocytes only to learn a hepatocyte-specific latent space,
then sweep Leiden resolutions for subclustering.

Input:
    Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad (657K hepatocytes)

Outputs (to results_gpu_v2/hepatocyte_subtypes/):
    hepatocyte_atlas.h5ad          — scVI latent + UMAP + all Leiden resolutions
    scvi_hepatocyte_model/         — saved scVI model weights
    leiden_resolution_scores.csv   — silhouette score per resolution

Environment: rapids_singlecell (GPU)
SBATCH: GPU partition, 1 GPU, 64GB RAM, 8 CPUs, --time=06:00:00
"""

import gc
import logging
import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
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
OUT_DIR = os.path.join(RESULTS, "hepatocyte_subtypes")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# GPU init
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.join(SC_DIR, "scripts"))
try:
    from gpu_utils import init_gpu, get_processor, to_gpu, from_gpu
    USE_GPU = init_gpu(pool_fraction=0.3)
except Exception:
    USE_GPU = False
    log.info("No GPU — running CPU-only")

# --- Seed pinning (added 2026-04-22 per T0.8) -------------------------------
import random
os.environ['PYTHONHASHSEED'] = '42'
random.seed(42)
np.random.seed(42)
import torch
torch.manual_seed(42)
torch.cuda.manual_seed_all(42)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
try:
    import scvi
    scvi.settings.seed = 42
except Exception:
    pass
# ---------------------------------------------------------------------------

import scanpy as sc
import anndata as ad
from scipy import sparse

if USE_GPU:
    pp, tl = get_processor(True)
else:
    pp, tl = sc.pp, sc.tl

# ---------------------------------------------------------------------------
# Fix h5ad categorical encoding (reused from 300)
# ---------------------------------------------------------------------------
def read_h5ad_fixed(h5ad_path):
    """Read h5ad with categorical 'ordered' attribute fix."""
    import h5py as _h5py
    try:
        return sc.read_h5ad(h5ad_path)
    except KeyError:
        pass
    import shutil
    tmp_path = os.path.join(OUT_DIR, ".atlas_patched.h5ad")
    if not os.path.exists(tmp_path):
        log.info("Patching h5ad categoricals (one-time copy)...")
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
    return sc.read_h5ad(tmp_path)


def main():
    log.info("=" * 60)
    log.info("308: Hepatocyte Subclustering — scVI Training")
    log.info("=" * 60)

    # ── 1. Load atlas and subset to hepatocytes ──────────────────────────
    log.info("Loading atlas from %s ...", H5AD)
    adata_full = read_h5ad_fixed(H5AD)
    log.info("Full atlas: %s", adata_full.shape)

    mask = adata_full.obs["cell_type"] == "Hepatocytes"
    log.info("Hepatocytes: %d cells", mask.sum())

    # Use raw counts for scVI (raw layer stores pre-normalization counts)
    # Must subset adata_full first (subsets both X and obs), then extract raw
    adata_hep = adata_full[mask].copy()
    del adata_full
    gc.collect()
    log.info("Subset to hepatocytes: %s", adata_hep.shape)

    if adata_hep.raw is not None:
        log.info("Using .raw layer for counts")
        raw_adata = adata_hep.raw.to_adata()
        # raw_adata has raw X + raw var but inherits obs from adata_hep
        adata = raw_adata
    else:
        log.info("No .raw layer — using X directly")
        adata = adata_hep

    del adata_hep

    gc.collect()
    log.info("Hepatocyte adata: %s", adata.shape)

    # Ensure counts are integers (scVI requires count data)
    if sparse.issparse(adata.X):
        x_min = adata.X.data.min() if adata.X.nnz > 0 else 0
        x_max = adata.X.data.max() if adata.X.nnz > 0 else 0
    else:
        x_min, x_max = adata.X.min(), adata.X.max()
    log.info("X range: [%.2f, %.2f]", x_min, x_max)

    # If X looks log-normalized (max < 20), we need raw counts
    # The .raw layer should have raw counts; verify by checking for integers
    if x_max < 20:
        log.warning("X appears log-normalized (max=%.2f). scVI needs raw counts.", x_max)
        log.warning("Attempting to undo log1p normalization...")
        # This is a fallback; .raw should have given us counts above
        if sparse.issparse(adata.X):
            adata.X = adata.X.expm1()
        else:
            adata.X = np.expm1(adata.X)
        log.info("After expm1 — X range: [%.2f, %.2f]",
                 adata.X.data.min() if sparse.issparse(adata.X) else adata.X.min(),
                 adata.X.data.max() if sparse.issparse(adata.X) else adata.X.max())

    # ── 2. Select HVGs (batch-aware) ────────────────────────────────────
    log.info("Selecting HVGs (batch_key='dataset', n_top=3000)...")

    # Store raw counts in a layer for scVI
    adata.layers["counts"] = adata.X.copy()

    # Normalize + log1p for HVG selection only
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Use 'seurat' flavor (log-normalized data) — 'seurat_v3' loess fit fails
    # on tiny batches (GSE136103: 192 cells, GSE189600: 175 cells)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=3000,
        batch_key="dataset",
        flavor="seurat",
    )

    # Exclude MT and ribosomal genes from HVGs
    mt_mask = adata.var_names.str.startswith("MT-")
    ribo_mask = adata.var_names.str.startswith("RPS") | adata.var_names.str.startswith("RPL")
    exclude = mt_mask | ribo_mask
    n_excluded = (adata.var["highly_variable"] & exclude).sum()
    adata.var.loc[exclude, "highly_variable"] = False
    log.info("Excluded %d MT/ribosomal genes from HVGs", n_excluded)
    n_hvg = adata.var["highly_variable"].sum()
    log.info("Final HVGs: %d", n_hvg)

    # Subset to HVGs — scVI does NOT auto-subset
    adata = adata[:, adata.var["highly_variable"]].copy()
    log.info("Subsetted to HVGs: %s", adata.shape)

    # ── 3. Train scVI ────────────────────────────────────────────────────
    log.info("Setting up scVI...")
    import scvi

    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key="dataset",
    )

    model = scvi.model.SCVI(
        adata,
        n_latent=20,
        n_layers=2,
        gene_likelihood="nb",
    )

    log.info("Training scVI (max 200 epochs)...")
    model.train(
        max_epochs=200,
        early_stopping=True,
        early_stopping_patience=10,
        train_size=0.9,
        batch_size=256,
    )

    # Save model
    model_dir = os.path.join(OUT_DIR, "scvi_hepatocyte_model")
    model.save(model_dir, overwrite=True)
    log.info("Saved scVI model to %s", model_dir)

    # ── 4. Get latent representation ─────────────────────────────────────
    log.info("Computing latent representation...")
    latent = model.get_latent_representation()
    adata.obsm["X_scVI"] = latent
    log.info("Latent shape: %s", latent.shape)

    # ── 5. Neighbors + UMAP ──────────────────────────────────────────────
    log.info("Computing neighbors (k=30) on scVI latent...")
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=30)

    log.info("Computing UMAP...")
    sc.tl.umap(adata)
    log.info("UMAP computed")

    # ── 6. Leiden sweep ──────────────────────────────────────────────────
    resolutions = [0.2, 0.3, 0.5, 0.7, 1.0]
    log.info("Sweeping Leiden resolutions: %s", resolutions)

    from sklearn.metrics import silhouette_score

    sil_scores = []
    for res in resolutions:
        key = f"leiden_{res}"
        sc.tl.leiden(adata, resolution=res, key_added=key)
        n_clusters = adata.obs[key].nunique()
        log.info("  res=%.1f → %d clusters", res, n_clusters)

        # Silhouette score on scVI latent (subsample for speed if > 50K cells)
        n_samp = min(50000, len(adata))
        np.random.seed(42)
        idx = np.random.choice(len(adata), n_samp, replace=False)
        sil = silhouette_score(
            latent[idx],
            adata.obs[key].values[idx],
            metric="euclidean",
            sample_size=None,
        )
        sil_scores.append({"resolution": res, "n_clusters": n_clusters, "silhouette": sil})
        log.info("  silhouette=%.4f", sil)

    sil_df = pd.DataFrame(sil_scores)
    sil_path = os.path.join(OUT_DIR, "leiden_resolution_scores.csv")
    sil_df.to_csv(sil_path, index=False)
    log.info("Saved resolution scores to %s", sil_path)

    # ── 7. Save h5ad ─────────────────────────────────────────────────────
    out_path = os.path.join(OUT_DIR, "hepatocyte_atlas.h5ad")
    log.info("Saving hepatocyte atlas to %s ...", out_path)
    adata.write_h5ad(out_path)
    log.info("Done. Shape: %s", adata.shape)


if __name__ == "__main__":
    main()
