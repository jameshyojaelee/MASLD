#!/usr/bin/env python3
"""
313: Compartment Subclustering via scVI (generic)

Clone of 308_hepatocyte_subcluster_scvi.py, generalized to subcluster ANY
pre-split compartment h5ad (endothelial, T cells, ...). Trains scVI on the
compartment to learn a compartment-specific latent space, then sweeps Leiden
resolutions.

KEY DIFFERENCE vs 308: the OUTPUT h5ad retains the FULL gene set (~37k genes),
not just the 3k HVGs. scVI is trained on an HVG sub-AnnData, but the latent
representation is written back onto the full-gene object so downstream receptor
queries (GLP1R, GIPR, GCGR, DPP4, GLP2R) and marker scoring have every gene.
The 308 template subsets to HVGs before writing, which would silently drop the
receptors we care about.

Input:
    a compartment h5ad with raw integer counts in .X (and/or layers['counts'])
    and obs columns: dataset (batch), condition, sample.

Outputs (to --outdir):
    {tag}_subcluster.h5ad          — scVI latent + UMAP + all Leiden resolutions,
                                      FULL gene set, log1p-CP10k in .X, raw in
                                      layers['counts']
    {tag}_scvi_model/              — saved scVI model weights
    {tag}_leiden_resolution_scores.csv — silhouette + n_clusters per resolution

Environment: rapids_singlecell (GPU)
SBATCH: --partition=gpu --gres=gpu:1 --mem=200G --cpus-per-task=16 --time=48:00:00

Usage:
    python 313_endothelial_subcluster_scvi.py \
        --input  .../atlas_cnmf_endothelial_cells.h5ad \
        --outdir .../results_gpu_v2/endothelial_subtypes \
        --tag    endothelial
"""

import argparse
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

# ---------------------------------------------------------------------------
# GPU init
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.join(SC_DIR, "scripts"))
try:
    from gpu_utils import init_gpu, get_processor
    USE_GPU = init_gpu(pool_fraction=0.3)
except Exception:
    USE_GPU = False
    log.info("No GPU — running CPU-only")

# --- Seed pinning (mirror 308 / T0.8) --------------------------------------
import random
os.environ["PYTHONHASHSEED"] = "42"
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
from scipy import sparse


def read_h5ad_fixed(h5ad_path, out_dir, tag):
    """Read h5ad with categorical 'ordered' attribute fix (reused from 308/300)."""
    import h5py as _h5py
    try:
        return sc.read_h5ad(h5ad_path)
    except KeyError:
        pass
    import shutil
    tmp_path = os.path.join(out_dir, f".{tag}_patched.h5ad")
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="compartment h5ad (raw counts)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tag", required=True, help="e.g. endothelial | tcell")
    ap.add_argument("--batch-key", default="dataset")
    ap.add_argument("--n-hvg", type=int, default=3000)
    ap.add_argument("--n-latent", type=int, default=20)
    ap.add_argument("--n-neighbors", type=int, default=30)
    ap.add_argument("--max-epochs", type=int, default=200)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    log.info("=" * 60)
    log.info("313: %s subclustering — scVI training", args.tag)
    log.info("=" * 60)

    # ── 1. Load compartment ──────────────────────────────────────────────
    log.info("Loading %s ...", args.input)
    adata = read_h5ad_fixed(args.input, args.outdir, args.tag)
    log.info("Compartment: %s (%d cells x %d genes)", args.tag, adata.n_obs, adata.n_vars)

    if args.batch_key not in adata.obs.columns:
        raise SystemExit(f"batch_key '{args.batch_key}' not in obs: {list(adata.obs.columns)}")

    # ── 2. Establish raw counts + log-normalized X ───────────────────────
    # Compartment h5ads differ: endothelial carries raw integer counts in .X +
    # layers['counts']; T-cell carries LOG-NORMALIZED .X (no counts layer) with
    # raw integer counts in .raw.X. scVI REQUIRES counts, so source them in
    # priority order: layers['counts'] -> .raw.X (gene-aligned) -> .X.
    # Canonicalize afterwards: layers['counts'] = raw (int), X = log1p-CP10k.
    if "counts" in adata.layers:
        log.info("Raw-count source: layers['counts']")
        raw = adata.layers["counts"]
    elif adata.raw is not None:
        log.info("Raw-count source: adata.raw.X (%d genes) — aligning to var_names",
                 adata.raw.n_vars)
        raw_ad = adata.raw.to_adata()
        missing = (~adata.var_names.isin(raw_ad.var_names)).sum()
        if missing:
            raise SystemExit(f"{missing} adata.var_names absent from .raw — cannot align counts")
        raw = raw_ad[:, adata.var_names].X  # reorder raw to match adata gene order
    else:
        log.info("Raw-count source: adata.X (no counts layer / no .raw)")
        raw = adata.X
    if not sparse.issparse(raw):
        raw = sparse.csr_matrix(raw)
    raw = raw.tocsr()
    raw.eliminate_zeros()

    xmax = raw.data.max() if raw.nnz else 0
    log.info("Raw counts range: [%.1f, %.1f], integer=%s",
             raw.data.min() if raw.nnz else 0, xmax,
             bool(np.allclose(raw.data[:100000], np.round(raw.data[:100000]))))
    if xmax < 20:
        log.warning("Raw counts max=%.1f looks log-normalized; scVI needs counts!", xmax)

    adata.layers["counts"] = raw.copy()
    # Build log1p-CP10k into .X (used for HVG + downstream marker/receptor scoring).
    adata.X = raw.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # ── 3. HVG selection (batch-aware) on a working copy ─────────────────
    log.info("Selecting HVGs (batch_key=%s, n_top=%d, flavor=seurat)...",
             args.batch_key, args.n_hvg)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=args.n_hvg,
        batch_key=args.batch_key,
        flavor="seurat",
    )
    # Exclude MT / ribosomal from HVGs
    mt = adata.var_names.str.startswith("MT-")
    ribo = adata.var_names.str.startswith("RPS") | adata.var_names.str.startswith("RPL")
    exclude = mt | ribo
    n_excl = int((adata.var["highly_variable"] & exclude).sum())
    adata.var.loc[exclude, "highly_variable"] = False
    n_hvg = int(adata.var["highly_variable"].sum())
    log.info("Excluded %d MT/ribo from HVGs; final HVGs=%d", n_excl, n_hvg)

    # ── 4. Train scVI on an HVG sub-AnnData (keep full adata intact) ──────
    import scvi
    model_adata = adata[:, adata.var["highly_variable"]].copy()
    log.info("Model AnnData (HVG subset): %s", model_adata.shape)

    scvi.model.SCVI.setup_anndata(
        model_adata,
        layer="counts",
        batch_key=args.batch_key,
    )
    model = scvi.model.SCVI(
        model_adata,
        n_latent=args.n_latent,
        n_layers=2,
        gene_likelihood="nb",
    )
    log.info("Training scVI (max %d epochs, early stop)...", args.max_epochs)
    model.train(
        max_epochs=args.max_epochs,
        early_stopping=True,
        early_stopping_patience=10,
        train_size=0.9,
        batch_size=256,
    )

    model_dir = os.path.join(args.outdir, f"{args.tag}_scvi_model")
    model.save(model_dir, overwrite=True)
    log.info("Saved scVI model to %s", model_dir)

    # ── 5. Latent → write back onto FULL-gene adata ──────────────────────
    latent = model.get_latent_representation()
    adata.obsm["X_scVI"] = latent
    log.info("Latent shape: %s (written onto full-gene adata: %s)", latent.shape, adata.shape)
    del model_adata
    gc.collect()

    # ── 6. Neighbors + UMAP (on scVI latent; CPU scanpy like 308) ────────
    log.info("Neighbors (k=%d) on X_scVI ...", args.n_neighbors)
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=args.n_neighbors)
    log.info("UMAP ...")
    sc.tl.umap(adata)

    # ── 7. Leiden sweep + silhouette ─────────────────────────────────────
    from sklearn.metrics import silhouette_score
    resolutions = [0.2, 0.3, 0.5, 0.7, 1.0]
    log.info("Leiden sweep: %s", resolutions)
    rows = []
    for res in resolutions:
        key = f"leiden_{res}"
        sc.tl.leiden(adata, resolution=res, key_added=key)
        n_clusters = int(adata.obs[key].nunique())
        n_samp = min(50000, adata.n_obs)
        np.random.seed(42)
        idx = np.random.choice(adata.n_obs, n_samp, replace=False)
        sil = float(silhouette_score(latent[idx], adata.obs[key].values[idx], metric="euclidean"))
        rows.append({"resolution": res, "n_clusters": n_clusters, "silhouette": sil})
        log.info("  res=%.1f -> %d clusters, silhouette=%.4f", res, n_clusters, sil)

    sil_df = pd.DataFrame(rows)
    sil_path = os.path.join(args.outdir, f"{args.tag}_leiden_resolution_scores.csv")
    sil_df.to_csv(sil_path, index=False)
    log.info("Saved resolution scores to %s", sil_path)

    # ── 8. Save FULL-gene h5ad ───────────────────────────────────────────
    out_path = os.path.join(args.outdir, f"{args.tag}_subcluster.h5ad")
    log.info("Saving %s ...", out_path)
    adata.write_h5ad(out_path)
    log.info("Done. Shape: %s", adata.shape)


if __name__ == "__main__":
    main()
