#!/usr/bin/env python
"""
02_train_scvi_v2.py — Phase 0.5 v2 atlas: scVI training

Trains scVI on the harmonized v2 hepatocyte atlas (Agent 2 output) with the
same hyperparameters as v1 Script 308 to enable head-to-head v1↔v2 comparison.

Inputs
------
- Analysis/SingleCell/results_gpu_v2_phase05/atlas/hepatocyte_atlas_v2.h5ad
    Agent 2's harmonized atlas. Expected to have:
      .layers["counts"] with raw integer counts
      .obs["dataset"]   batch key (per-dataset)
      .obs["sample"]    donor key (per-donor)
    The atlas may or may not be HVG-subset; we re-select HVGs here (seurat_v3,
    batch-aware, ~3000 genes) so the v2 HVG list is reproducible from a single
    script.

Outputs (all under Analysis/SingleCell/results_gpu_v2_phase05/)
---------------------------------------------------------------
- scvi/scvi_hepatocyte_model_v2/        — saved scVI model (model.pt + config)
- scvi/canonical_hvg_genes_v2.tsv       — frozen HVG list (per Phase 4 spec)
- atlas/hepatocyte_atlas_v2_annotated.h5ad — atlas with X_scVI in .obsm
- scvi/scvi_training_log.json           — hyperparameters + final epoch + ELBO

Hyperparameters (match v1 Script 308 for head-to-head comparability)
--------------------------------------------------------------------
- batch_key='dataset'
- n_latent=20, n_layers=2, n_hidden=128
- gene_likelihood='nb'
- dropout_rate=0.1
- max_epochs=150 with early stopping (patience=10)
- train_size=0.9, batch_size=512

Environment: .agent/scanpy_env (Python 3.10, torch 2.3.1+cu121, scvi-tools 1.3.3)
SLURM: gpu partition, L40S (exclude Blackwell ne1dg7-[001-010]), 1 GPU, 64GB, 8 CPUs, 24h
"""
from __future__ import annotations

import gc
import json
import logging
import os
import sys
import time
import warnings
from pathlib import Path

# scvi-tools network / telemetry off
os.environ.setdefault("WANDB_MODE", "disabled")
os.environ.setdefault("WANDB_DISABLED", "true")
os.environ.setdefault("SCVI_NO_TELEMETRY", "true")

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

# Seed pinning (matches v1 Script 308)
import random
os.environ["PYTHONHASHSEED"] = "42"
random.seed(42)
np.random.seed(42)
import torch
torch.manual_seed(42)
torch.cuda.manual_seed_all(42)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

import anndata as ad
import scanpy as sc
import scvi
from scipy import sparse

scvi.settings.seed = 42

# ── Paths ───────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
PHASE05 = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2_phase05"

ATLAS_IN = PHASE05 / "atlas/hepatocyte_atlas_v2.h5ad"
ATLAS_OUT = PHASE05 / "atlas/hepatocyte_atlas_v2_annotated.h5ad"
SCVI_DIR = PHASE05 / "scvi"
MODEL_OUT = SCVI_DIR / "scvi_hepatocyte_model_v2"
HVG_TSV = SCVI_DIR / "canonical_hvg_genes_v2.tsv"
TRAIN_LOG = SCVI_DIR / "scvi_training_log.json"

# Hyperparameters
SCVI_N_LATENT = 20
SCVI_N_LAYERS = 2
SCVI_N_HIDDEN = 128
SCVI_GENE_LIKELIHOOD = "nb"
SCVI_DROPOUT = 0.1
SCVI_MAX_EPOCHS = 150
SCVI_PATIENCE = 10
SCVI_TRAIN_SIZE = 0.9
SCVI_BATCH_SIZE = 512
HVG_N_TOP = 3000
HVG_FLAVOR = "seurat_v3"


def _register_null_io_fallback():
    """anndata>=0.10 may write encoding_type='null'; older readers crash."""
    try:
        import h5py
        from anndata._io.specs.registry import _REGISTRY, IOSpec
        _REGISTRY.register_read(h5py.Dataset, IOSpec("null", "0.1.0"))(
            lambda elem, _reader=None: None)
        _REGISTRY.register_read(h5py.Group, IOSpec("null", "0.1.0"))(
            lambda elem, _reader=None: None)
    except Exception as e:
        log.warning(f"failed to register null IO fallback ({e})")


def load_atlas() -> ad.AnnData:
    _register_null_io_fallback()
    log.info(f"loading {ATLAS_IN}")
    if not ATLAS_IN.exists():
        log.error(f"FATAL: atlas missing at {ATLAS_IN}")
        sys.exit(1)
    a = ad.read_h5ad(ATLAS_IN)
    log.info(f"  shape: {a.shape}")
    log.info(f"  obs cols: {list(a.obs.columns)[:20]}")
    log.info(f"  obsm: {list(a.obsm.keys())}")
    log.info(f"  layers: {list(a.layers.keys())}")
    return a


def ensure_counts_layer(adata: ad.AnnData) -> ad.AnnData:
    """Make sure adata.layers['counts'] holds raw integer counts."""
    if "counts" in adata.layers:
        log.info("'counts' layer present")
        return adata

    # Heuristic: if .raw is available, use raw counts
    if adata.raw is not None:
        log.info("no 'counts' layer; pulling from adata.raw")
        raw_adata = adata.raw.to_adata()
        # Keep obs/obsm from adata, X from raw
        adata.layers["counts"] = raw_adata.X
        return adata

    # Fallback: if X looks like raw counts (integers), use directly
    if sparse.issparse(adata.X):
        x_max = adata.X.data.max() if adata.X.nnz > 0 else 0
        x_min = adata.X.data.min() if adata.X.nnz > 0 else 0
    else:
        x_max = float(adata.X.max())
        x_min = float(adata.X.min())
    log.info(f"X range: [{x_min:.2f}, {x_max:.2f}]")

    if x_max > 20 and x_min >= 0:
        log.info("X looks like raw counts; copying to 'counts' layer")
        adata.layers["counts"] = adata.X.copy()
    else:
        log.error(
            "FATAL: cannot find raw counts. atlas v2 must provide either "
            "layers['counts'], adata.raw, or raw integer X."
        )
        sys.exit(1)
    return adata


def select_hvgs(adata: ad.AnnData) -> ad.AnnData:
    """Run seurat_v3 HVG (batch-aware) on the raw-count layer.

    Per Phase 1 spec: batch_key='dataset', flavor='seurat_v3', n_top=3000,
    MT/ribo excluded from HVGs.
    """
    log.info(f"selecting HVGs: flavor={HVG_FLAVOR}, n_top={HVG_N_TOP}, batch_key='dataset'")

    # Pre-filter near-zero-variance / ultra-rare genes that destabilize the
    # seurat_v3 LOESS fit. Light filter (min_cells=10) on the global matrix
    # — per-batch LOESS still tolerates genes that are 0-variance in 1-2 of 7
    # datasets because seurat_v3 takes the median rank across batches.
    n_genes_before = adata.shape[1]
    # filter_genes uses .X by default — point it at the counts layer
    cell_counts_per_gene = np.asarray((adata.layers["counts"] > 0).sum(axis=0)).ravel()
    keep = cell_counts_per_gene >= 10
    adata = adata[:, keep].copy()
    log.info(
        f"  min_cells=10 filter (on counts layer): {n_genes_before} -> {adata.shape[1]} genes"
    )

    # seurat_v3 wants raw counts, not log-normalized. Fall back to 'seurat'
    # flavor if LOESS still produces a near-singular fit.
    try:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=HVG_N_TOP,
            batch_key="dataset",
            flavor=HVG_FLAVOR,
            layer="counts",
        )
    except ValueError as e:
        log.warning(f"seurat_v3 HVG failed ({e}); falling back to 'seurat' on log1p")
        # 'seurat' flavor needs log1p data in .X — make a temporary copy
        tmp = adata.copy()
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
        sc.pp.highly_variable_genes(
            tmp,
            n_top_genes=HVG_N_TOP,
            batch_key="dataset",
            flavor="seurat",
        )
        adata.var["highly_variable"] = tmp.var["highly_variable"].values
        if "highly_variable_rank" in tmp.var.columns:
            adata.var["highly_variable_rank"] = tmp.var["highly_variable_rank"].values
        del tmp

    # Exclude MT/ribo
    mt = adata.var_names.str.startswith("MT-")
    ribo = adata.var_names.str.startswith("RPS") | adata.var_names.str.startswith("RPL")
    exclude = mt | ribo
    n_excluded = int((adata.var["highly_variable"] & exclude).sum())
    adata.var.loc[exclude, "highly_variable"] = False
    n_hvg = int(adata.var["highly_variable"].sum())
    log.info(f"  excluded {n_excluded} MT/ribo from HVGs; final HVG count: {n_hvg}")

    # Save canonical HVG list BEFORE subsetting
    SCVI_DIR.mkdir(parents=True, exist_ok=True)
    hvg_list = pd.DataFrame({
        "gene_symbol": adata.var_names[adata.var["highly_variable"]].astype(str),
        "highly_variable_rank": np.argsort(
            -adata.var.loc[adata.var["highly_variable"], "highly_variable_rank"].fillna(0).values
        ).argsort() if "highly_variable_rank" in adata.var.columns else range(n_hvg),
    })
    hvg_list.to_csv(HVG_TSV, sep="\t", index=False)
    log.info(f"  saved HVG list -> {HVG_TSV}")

    # Subset to HVGs for scVI
    adata = adata[:, adata.var["highly_variable"]].copy()
    log.info(f"subset to HVGs: {adata.shape}")
    return adata


def train_scvi(adata: ad.AnnData) -> "scvi.model.SCVI":
    log.info("setting up scvi-tools AnnData")
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key="dataset",
    )
    log.info("building scVI model")
    model = scvi.model.SCVI(
        adata,
        n_latent=SCVI_N_LATENT,
        n_layers=SCVI_N_LAYERS,
        n_hidden=SCVI_N_HIDDEN,
        gene_likelihood=SCVI_GENE_LIKELIHOOD,
        dropout_rate=SCVI_DROPOUT,
    )
    log.info(f"training scVI (max_epochs={SCVI_MAX_EPOCHS}, patience={SCVI_PATIENCE})")
    t0 = time.time()
    model.train(
        max_epochs=SCVI_MAX_EPOCHS,
        early_stopping=True,
        early_stopping_patience=SCVI_PATIENCE,
        train_size=SCVI_TRAIN_SIZE,
        batch_size=SCVI_BATCH_SIZE,
        enable_progress_bar=False,
    )
    elapsed = time.time() - t0
    log.info(f"  training done in {elapsed/60:.1f} min")
    return model


def main():
    log.info("=" * 70)
    log.info("02_train_scvi_v2: scVI v2 training on harmonized atlas")
    log.info("=" * 70)

    log.info(f"GPU available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        log.info(f"  device: {torch.cuda.get_device_name(0)}")
        log.info(f"  torch version: {torch.__version__}")
    else:
        log.warning("NO GPU detected — training will be impractically slow")

    MODEL_OUT.mkdir(parents=True, exist_ok=True)
    ATLAS_OUT.parent.mkdir(parents=True, exist_ok=True)

    adata = load_atlas()

    # Validate required obs columns
    for col in ("dataset", "sample"):
        if col not in adata.obs.columns:
            log.error(f"FATAL: required obs column '{col}' missing from atlas v2")
            sys.exit(1)
    log.info(
        f"  datasets: {adata.obs['dataset'].nunique()}; "
        f"samples (donors): {adata.obs['sample'].nunique()}"
    )

    adata = ensure_counts_layer(adata)

    # If atlas v2 is already HVG-subset (var has 'highly_variable' all True and
    # is small), skip HVG selection. Otherwise compute HVGs.
    needs_hvg = True
    if "highly_variable" in adata.var.columns and adata.n_vars <= 5000:
        already = int(adata.var["highly_variable"].sum())
        if already > 0 and already >= adata.n_vars * 0.9:
            log.info(
                f"atlas appears pre-HVG-subset ({adata.n_vars} genes, "
                f"{already} marked HVG); skipping HVG re-selection"
            )
            needs_hvg = False
            # Persist the existing HVG list anyway for reproducibility
            SCVI_DIR.mkdir(parents=True, exist_ok=True)
            pd.DataFrame({
                "gene_symbol": adata.var_names.astype(str),
            }).to_csv(HVG_TSV, sep="\t", index=False)
            log.info(f"  persisted existing HVG list -> {HVG_TSV}")
    if needs_hvg:
        adata = select_hvgs(adata)

    log.info(f"final adata for scVI: {adata.shape}")
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

    # Train scVI
    model = train_scvi(adata)

    # Save model
    log.info(f"saving scVI model to {MODEL_OUT}")
    model.save(str(MODEL_OUT), overwrite=True)

    # Inject latent into atlas obsm
    log.info("computing latent representation")
    latent = model.get_latent_representation()
    adata.obsm["X_scVI"] = latent
    log.info(f"  latent shape: {latent.shape}")

    # Re-load the FULL atlas (not HVG-subset) and stamp X_scVI on it so
    # downstream scANVI / UMAP / DE work on the full gene set.
    log.info("re-loading full atlas to stamp X_scVI obsm")
    full = load_atlas()
    # Sanity: cell order must match. If not, sort latent by adata.obs_names.
    if list(full.obs_names) != list(adata.obs_names):
        log.info("  cell order differs; reindexing latent to full atlas order")
        order = pd.Index(adata.obs_names).get_indexer(full.obs_names)
        if (order == -1).any():
            log.error("FATAL: subset cells missing in full atlas — cannot align")
            sys.exit(1)
        latent_full = latent[order]
    else:
        latent_full = latent
    full.obsm["X_scVI"] = latent_full

    # Persist `highly_variable` so downstream knows which genes were used
    full.var["highly_variable_v2"] = full.var_names.isin(adata.var_names)

    log.info(f"writing annotated atlas -> {ATLAS_OUT}")
    full.write_h5ad(ATLAS_OUT, compression="gzip")
    log.info(f"  shape: {full.shape}; obsm: {list(full.obsm.keys())}")

    # Training log
    train_record = {
        "script": "02_train_scvi_v2.py",
        "atlas_in": str(ATLAS_IN),
        "atlas_out": str(ATLAS_OUT),
        "model_dir": str(MODEL_OUT),
        "hvg_tsv": str(HVG_TSV),
        "n_cells": int(adata.n_obs),
        "n_hvgs": int(adata.n_vars),
        "n_datasets": int(adata.obs["dataset"].nunique()),
        "n_donors": int(adata.obs["sample"].nunique()),
        "hyperparams": {
            "n_latent": SCVI_N_LATENT,
            "n_layers": SCVI_N_LAYERS,
            "n_hidden": SCVI_N_HIDDEN,
            "gene_likelihood": SCVI_GENE_LIKELIHOOD,
            "dropout_rate": SCVI_DROPOUT,
            "max_epochs": SCVI_MAX_EPOCHS,
            "patience": SCVI_PATIENCE,
            "train_size": SCVI_TRAIN_SIZE,
            "batch_size": SCVI_BATCH_SIZE,
            "hvg_flavor": HVG_FLAVOR,
            "hvg_n_top": HVG_N_TOP,
        },
        "scvi_tools_version": scvi.__version__,
        "torch_version": torch.__version__,
        "cuda_available": bool(torch.cuda.is_available()),
        "device_name": (
            torch.cuda.get_device_name(0) if torch.cuda.is_available() else "cpu"
        ),
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    with open(TRAIN_LOG, "w") as f:
        json.dump(train_record, f, indent=2)
    log.info(f"wrote training log -> {TRAIN_LOG}")
    log.info("DONE")


if __name__ == "__main__":
    main()
