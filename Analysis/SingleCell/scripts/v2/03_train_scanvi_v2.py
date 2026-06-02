#!/usr/bin/env python
"""
03_train_scanvi_v2.py — Phase 0.5 v2 atlas: scANVI fine-tuning (two variants)

Trains scANVI on top of the persisted scVI v2 base model with two label-broadcast
strategies, runs 5-fold CV across Andrews donors, and writes per-cell + per-donor
predictions. Designed to mirror v1 scripts 343h (vanilla) and 343n (augmented)
so the v2 → v1 comparison is apples-to-apples.

Variants (selected via --variant)
---------------------------------
- vanilla   : broadcast Andrews documented F-stage labels only (~75K labelled cells)
- augmented : broadcast Andrews documented + clean-healthy donors → F0 + Cirrhosis → F4
              (~186K labelled cells; matches Script 343n)

Inputs
------
- Analysis/SingleCell/results_gpu_v2_phase05/atlas/hepatocyte_atlas_v2_annotated.h5ad
    Agent 3 step-1 output (this script's step-1 sibling); has X_scVI in obsm.
- Analysis/SingleCell/results_gpu_v2_phase05/scvi/scvi_hepatocyte_model_v2/model.pt
    Persisted scVI v2 base model.
- Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_documented.tsv
    58 Andrews donors with documented Kleiner F-stage (read-only, shared with v1).
- Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv
    Donor coarse stage labels (read-only).
- Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/dubious_healthy_donors.tsv
    Donors flagged as not-clean-healthy (read-only).

Outputs (under Analysis/SingleCell/results_gpu_v2_phase05/scvi/)
-----------------------------------------------------------------
- scanvi_<variant>_v2_model/scanvi_full/ — persisted scANVI model
- donor_fstage_scanvi_<variant>_v2_predicted.tsv
- percell_fstage_scanvi_<variant>_v2.tsv.gz
- scanvi_<variant>_v2_method_metrics.tsv
- scanvi_<variant>_v2_external_validation_summary.txt
- scanvi_<variant>_v2_training_log.json

Hyperparameters
---------------
- n_samples_per_label = 200 (matches Script 343n)
- max_epochs = 75
- 5-fold StratifiedKFold over Andrews donors, seed=42

Environment: .agent/scanpy_env (Python 3.10, torch 2.3.1+cu121, scvi-tools 1.3.3)
SLURM: gpu partition, L40S (exclude Blackwell ne1dg7-[001-010]), 1 GPU, 64GB, 8 CPUs, 24h
"""
from __future__ import annotations

import argparse
import gc
import json
import logging
import os
import sys
import time
import warnings
from collections import Counter
from pathlib import Path

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
import scvi
from scipy import stats, sparse
from sklearn.metrics import (
    accuracy_score,
    cohen_kappa_score,
    confusion_matrix,
)
from sklearn.model_selection import StratifiedKFold

scvi.settings.seed = 42


# ── Paths ───────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
PHASE05 = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2_phase05"

ATLAS_IN = PHASE05 / "atlas/hepatocyte_atlas_v2_annotated.h5ad"
SCVI_MODEL_DIR = PHASE05 / "scvi/scvi_hepatocyte_model_v2"
SCVI_DIR = PHASE05 / "scvi"

# Read-only v1 metadata (donor F-stage docs, coarse stage, dubious-healthy)
V1_STAGE_DIR = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
FSTAGE_DOC = V1_STAGE_DIR / "donor_fstage_documented.tsv"
META_EXT = V1_STAGE_DIR / "donor_metadata_extended.tsv"
DUBIOUS = V1_STAGE_DIR / "dubious_healthy_donors.tsv"

# Hyperparameters
SCVI_N_LATENT = 20
SCVI_N_LAYERS = 2
SCVI_N_HIDDEN = 128
SCVI_GENE_LIKELIHOOD = "nb"
SCVI_DROPOUT = 0.1

SCANVI_MAX_EPOCHS = 75
SCANVI_N_SAMPLES_PER_LABEL = 200  # matches Script 343n
CV_N_SPLITS = 5

UNLABELED = "Unknown"
F_STAGES = [0, 1, 2, 3, 4]
F_STAGE_STR = [str(f) for f in F_STAGES]


def _register_null_io_fallback():
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
        log.error(f"FATAL: annotated atlas missing at {ATLAS_IN}")
        sys.exit(1)
    a = ad.read_h5ad(ATLAS_IN)
    log.info(f"  shape: {a.shape}; obsm: {list(a.obsm.keys())}; layers: {list(a.layers.keys())}")
    if "counts" not in a.layers:
        log.error("FATAL: 'counts' layer missing on annotated atlas")
        sys.exit(1)
    # If atlas has highly_variable_v2 flag, subset to HVG genes (scVI was trained on those)
    if "highly_variable_v2" in a.var.columns:
        n_hvg = int(a.var["highly_variable_v2"].sum())
        log.info(f"  subsetting to {n_hvg} v2 HVG genes")
        a = a[:, a.var["highly_variable_v2"]].copy()
    return a


def broadcast_vanilla_labels(adata: ad.AnnData, fstage_df: pd.DataFrame) -> pd.Series:
    """Andrews documented donors only → cell label = F-stage 0..4; else Unknown."""
    fstage_df = fstage_df.dropna(subset=["F_stage_documented"]).copy()
    fstage_df["F_stage_documented"] = pd.to_numeric(
        fstage_df["F_stage_documented"], errors="coerce")
    fstage_df = fstage_df.dropna(subset=["F_stage_documented"])
    fstage_df["F_stage_documented"] = fstage_df["F_stage_documented"].astype(int)

    donor_to_f = dict(zip(fstage_df["sample"], fstage_df["F_stage_documented"]))
    samples = adata.obs["sample"].astype(str).values
    labels = np.array([
        str(int(donor_to_f[s])) if s in donor_to_f else UNLABELED for s in samples
    ], dtype=object)
    labels = pd.Series(labels, index=adata.obs_names)
    n_lab = int((labels != UNLABELED).sum())
    log.info(f"  vanilla broadcast: {n_lab:,}/{len(labels):,} cells labelled")
    log.info(f"    per-class counts: {dict(Counter(labels))}")
    return labels


def broadcast_augmented_labels(
    adata: ad.AnnData,
    fstage_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    dubious_samples: set,
) -> pd.Series:
    """Andrews + clean-Healthy → F0 + Cirrhosis → F4."""
    fstage_df = fstage_df.dropna(subset=["F_stage_documented"]).copy()
    fstage_df["F_stage_documented"] = pd.to_numeric(
        fstage_df["F_stage_documented"], errors="coerce")
    fstage_df = fstage_df.dropna(subset=["F_stage_documented"])
    fstage_df["F_stage_documented"] = fstage_df["F_stage_documented"].astype(int)
    docs = dict(zip(fstage_df["sample"], fstage_df["F_stage_documented"]))

    coarse = dict(zip(meta_df["sample"], meta_df["disease_stage_coarse"]))

    samples = adata.obs["sample"].astype(str).values
    labels = np.full(len(samples), UNLABELED, dtype=object)
    n_doc = n_chealthy = n_cirr = 0
    for i, s in enumerate(samples):
        if s in docs:
            labels[i] = str(int(docs[s]))
            n_doc += 1
        else:
            c = coarse.get(s, None)
            if c == "Healthy" and (s not in dubious_samples):
                labels[i] = "0"
                n_chealthy += 1
            elif c == "Cirrhosis":
                labels[i] = "4"
                n_cirr += 1
    labels = pd.Series(labels, index=adata.obs_names)
    n_lab = int((labels != UNLABELED).sum())
    log.info(f"  augmented broadcast: {n_lab:,}/{len(labels):,} cells labelled")
    log.info(f"    Andrews documented: {n_doc:,}")
    log.info(f"    clean-Healthy -> F0: {n_chealthy:,}")
    log.info(f"    Cirrhosis -> F4    : {n_cirr:,}")
    log.info(f"    per-class counts: {dict(Counter(labels))}")
    return labels


def predict_with_proba(scanvi_model, adata) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Return (argmax_int, P[N,5], class_labels_str) in F0..F4 order."""
    soft = scanvi_model.predict(adata, soft=True)
    if UNLABELED in soft.columns:
        soft = soft.drop(columns=[UNLABELED])
    for s in F_STAGE_STR:
        if s not in soft.columns:
            soft[s] = 0.0
    soft = soft[F_STAGE_STR]
    rowsum = soft.values.sum(axis=1, keepdims=True)
    rowsum[rowsum == 0] = 1.0
    proba = soft.values / rowsum
    argmax = proba.argmax(axis=1)
    return argmax, proba, F_STAGE_STR


def aggregate_donor(per_cell: pd.DataFrame, p_cols: list, argmax_col: str) -> pd.DataFrame:
    """Donor-level aggregation: argmax of mean posterior + mode of per-cell argmax."""
    grouped_mean = per_cell.groupby(["sample", "dataset"], observed=True)[p_cols].mean()
    grouped_mean.columns = [
        f"posterior_F{c.split('P_F')[-1]}_mean" for c in p_cols
    ]
    argmax_per_donor = grouped_mean.values.argmax(axis=1)
    mode_per_donor = (
        per_cell.groupby(["sample", "dataset"], observed=True)[argmax_col]
        .agg(lambda x: x.value_counts().idxmax())
    )
    n_hep = per_cell.groupby(["sample", "dataset"], observed=True).size().rename("n_cells")
    posterior_mean = (
        grouped_mean.values * np.array([0, 1, 2, 3, 4])
    ).sum(axis=1)
    out = grouped_mean.copy()
    out[f"{argmax_col}_argmax"] = argmax_per_donor
    out[f"{argmax_col}_mode"] = mode_per_donor.values
    out["posterior_mean"] = posterior_mean
    out["n_cells"] = n_hep.values
    return out.reset_index()


def fold_cv(
    adata: ad.AnnData,
    fstage_df: pd.DataFrame,
    full_cell_labels: pd.Series,
    scvi_model_main,
    labels_col_prefix: str,
    argmax_col: str,
    n_splits: int = 5,
    seed: int = 42,
) -> dict:
    """5-fold stratified CV across Andrews documented donors.

    Each fold relabels held-out donors' cells as 'Unknown' and re-fits the
    scANVI head from the frozen scVI base.
    """
    docs = fstage_df.dropna(subset=["F_stage_documented"]).copy()
    docs["F_stage_documented"] = docs["F_stage_documented"].astype(int)
    donors = docs["sample"].values
    y = docs["F_stage_documented"].values

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    fold_preds = []
    full_arr = full_cell_labels.values.copy()
    cell_samples = adata.obs["sample"].astype(str).values

    for fold_idx, (train_idx, test_idx) in enumerate(skf.split(donors, y)):
        held_out = set(donors[test_idx])
        log.info(f"FOLD {fold_idx+1}/{n_splits}: holding out {len(held_out)} Andrews donors")

        fold_labels = full_arr.copy()
        held_mask = np.isin(cell_samples, list(held_out))
        fold_labels[held_mask] = UNLABELED
        n_held = int(held_mask.sum())
        n_train_lab = int((fold_labels != UNLABELED).sum())
        log.info(f"  cells held out: {n_held:,}; training-labelled: {n_train_lab:,}")

        col = f"{labels_col_prefix}_fold{fold_idx}"
        adata.obs[col] = pd.Categorical(
            fold_labels, categories=F_STAGE_STR + [UNLABELED])

        scanvi = scvi.model.SCANVI.from_scvi_model(
            scvi_model_main,
            labels_key=col,
            unlabeled_category=UNLABELED,
            adata=adata,
        )
        scanvi.train(
            max_epochs=SCANVI_MAX_EPOCHS,
            n_samples_per_label=SCANVI_N_SAMPLES_PER_LABEL,
            check_val_every_n_epoch=10,
            enable_progress_bar=False,
        )
        a_held = adata[held_mask]
        argmax_int, proba, classes = predict_with_proba(scanvi, a_held)
        per_cell = pd.DataFrame({
            "sample": a_held.obs["sample"].astype(str).values,
            "dataset": a_held.obs["dataset"].astype(str).values,
            argmax_col: argmax_int,
        })
        for j, c in enumerate(classes):
            per_cell[f"P_F{c}"] = proba[:, j]
        p_cols = [f"P_F{c}" for c in classes]
        donor_agg = aggregate_donor(per_cell, p_cols, argmax_col)
        donor_agg["fold"] = fold_idx
        fold_preds.append(donor_agg)

        del scanvi
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    cv = pd.concat(fold_preds, ignore_index=True)
    cv = cv.merge(docs[["sample", "F_stage_documented"]], on="sample", how="left")
    qwk = cohen_kappa_score(
        cv["F_stage_documented"].astype(int),
        cv[f"{argmax_col}_argmax"].astype(int),
        weights="quadratic", labels=F_STAGES,
    )
    acc = accuracy_score(
        cv["F_stage_documented"].astype(int),
        cv[f"{argmax_col}_argmax"].astype(int),
    )
    cm = confusion_matrix(
        cv["F_stage_documented"].astype(int),
        cv[f"{argmax_col}_argmax"].astype(int),
        labels=F_STAGES,
    )
    log.info(f"{n_splits}-fold CV: QWK={qwk:.3f}, accuracy={acc:.3f}")
    log.info(f"confusion matrix (rows=true F-stage, cols=pred):\n{cm}")
    return dict(qwk=float(qwk), accuracy=float(acc), cm=cm.tolist(),
                cv=cv, n_splits=n_splits)


def external_validation(donor_preds: pd.DataFrame, meta: pd.DataFrame,
                        argmax_col: str) -> dict:
    df = donor_preds.merge(
        meta[["sample", "disease_stage_coarse", "disease_stage_numeric"]],
        on="sample", how="left",
    )
    df = df.dropna(subset=[f"{argmax_col}_argmax", "disease_stage_numeric"])
    rho_overall, p_overall = stats.spearmanr(
        df[f"{argmax_col}_argmax"], df["disease_stage_numeric"]
    )
    log.info(f"external Spearman rho({argmax_col}, disease_stage_numeric)"
             f" = {rho_overall:.3f} (n={len(df)})")
    per_ds = []
    for ds, g in df.groupby("dataset"):
        if g["disease_stage_numeric"].nunique() < 2:
            per_ds.append(dict(dataset=ds, rho=np.nan, p=np.nan, n=int(len(g))))
            continue
        r, p = stats.spearmanr(g[f"{argmax_col}_argmax"], g["disease_stage_numeric"])
        per_ds.append(dict(dataset=ds, rho=float(r), p=float(p), n=int(len(g))))
    per_ds_df = pd.DataFrame(per_ds)
    n_healthy_f4 = int(((df[f"{argmax_col}_argmax"] == 4) &
                        (df["disease_stage_coarse"] == "Healthy")).sum())
    cross = pd.crosstab(
        df[f"{argmax_col}_argmax"], df["disease_stage_coarse"]
    ).reindex(
        index=F_STAGES,
        columns=["Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"],
        fill_value=0,
    )
    log.info(f"Healthy donors predicted F4: {n_healthy_f4}")
    log.info(f"predicted F-stage x disease_stage_coarse:\n{cross}")
    return dict(
        rho_overall=float(rho_overall),
        p_overall=float(p_overall),
        per_dataset=per_ds_df,
        n_healthy_f4=n_healthy_f4,
        crosstab=cross,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--variant", required=True, choices=["vanilla", "augmented"],
        help="scANVI label-broadcast strategy",
    )
    ap.add_argument("--skip-cv", action="store_true",
                    help="Skip 5-fold CV (just train once and predict).")
    ap.add_argument("--cv-folds", type=int, default=CV_N_SPLITS)
    args = ap.parse_args()

    variant = args.variant
    log.info("=" * 70)
    log.info(f"03_train_scanvi_v2: scANVI v2 ({variant}) training")
    log.info("=" * 70)

    # variant-specific output paths
    MODEL_OUT = SCVI_DIR / f"scanvi_{variant}_v2_model"
    OUT_DONOR_TSV = SCVI_DIR / f"donor_fstage_scanvi_{variant}_v2_predicted.tsv"
    OUT_PERCELL = SCVI_DIR / f"percell_fstage_scanvi_{variant}_v2.tsv.gz"
    OUT_METRICS = SCVI_DIR / f"scanvi_{variant}_v2_method_metrics.tsv"
    OUT_SUMMARY = SCVI_DIR / f"scanvi_{variant}_v2_external_validation_summary.txt"
    OUT_LOG = SCVI_DIR / f"scanvi_{variant}_v2_training_log.json"
    MODEL_OUT.mkdir(parents=True, exist_ok=True)

    labels_col = f"fstage_{variant}_label"
    argmax_col = f"F_stage_scanvi_{variant}"

    log.info(f"GPU available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        log.info(f"  device: {torch.cuda.get_device_name(0)}")
    else:
        log.warning("NO GPU detected — training will be impractically slow")

    adata = load_atlas()

    # Validate required obs columns
    for col in ("dataset", "sample"):
        if col not in adata.obs.columns:
            log.error(f"FATAL: required obs column '{col}' missing")
            sys.exit(1)

    # ── Load broadcasted labels ─────────────────────────────────────────────
    log.info(f"loading donor F-stage metadata from v1 stage_trajectory dir")
    fstage = pd.read_csv(FSTAGE_DOC, sep="\t")
    fstage["F_stage_documented"] = pd.to_numeric(
        fstage["F_stage_documented"], errors="coerce")
    meta = pd.read_csv(META_EXT, sep="\t")

    if variant == "vanilla":
        cell_labels = broadcast_vanilla_labels(adata, fstage)
    else:
        dubious = pd.read_csv(DUBIOUS, sep="\t")
        dubious_samples = set(dubious["sample"].astype(str).values)
        log.info(f"  {len(dubious_samples)} dubious-healthy donors loaded")
        cell_labels = broadcast_augmented_labels(adata, fstage, meta, dubious_samples)

    n_labelled = int((cell_labels != UNLABELED).sum())
    adata.obs[labels_col] = pd.Categorical(
        cell_labels, categories=F_STAGE_STR + [UNLABELED])

    # ── Load persisted scVI base ───────────────────────────────────────────
    log.info(f"loading persisted scVI v2 from {SCVI_MODEL_DIR}")
    if not (SCVI_MODEL_DIR / "model.pt").exists():
        log.error(f"FATAL: scVI v2 model missing at {SCVI_MODEL_DIR / 'model.pt'}")
        sys.exit(1)
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="dataset")
    try:
        scvi_model = scvi.model.SCVI.load(str(SCVI_MODEL_DIR), adata=adata)
        log.info("  loaded scVI v2")
    except Exception as e:
        log.error(f"FATAL: could not load scVI v2 model ({e})")
        sys.exit(1)

    # ── 5-fold CV (across Andrews documented donors) ───────────────────────
    cv_result = None
    if not args.skip_cv:
        log.info("===== 5-fold CV =====")
        try:
            cv_result = fold_cv(
                adata, fstage, cell_labels, scvi_model,
                labels_col_prefix=f"fstage_{variant}",
                argmax_col=argmax_col,
                n_splits=args.cv_folds, seed=42,
            )
        except Exception as e:
            log.exception(f"CV FAILED: {e}")
            cv_result = None

    # ── Final scANVI: train on full broadcast labels, predict all cells ─────
    log.info("===== final scANVI on full broadcast =====")
    # Reset labels to the full broadcast
    adata.obs[labels_col] = pd.Categorical(
        cell_labels, categories=F_STAGE_STR + [UNLABELED])
    scanvi_full = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        labels_key=labels_col,
        unlabeled_category=UNLABELED,
        adata=adata,
    )
    scanvi_full.train(
        max_epochs=SCANVI_MAX_EPOCHS,
        n_samples_per_label=SCANVI_N_SAMPLES_PER_LABEL,
        check_val_every_n_epoch=10,
        enable_progress_bar=False,
    )
    scanvi_full.save(str(MODEL_OUT / "scanvi_full"), overwrite=True)
    log.info(f"  saved scANVI model -> {MODEL_OUT / 'scanvi_full'}")

    # ── Per-cell predictions for all hepatocytes ───────────────────────────
    log.info("predicting F-stage for all hepatocytes")
    argmax_int, proba, classes = predict_with_proba(scanvi_full, adata)
    per_cell = pd.DataFrame({
        "cell_id": adata.obs_names.values.astype(str),
        "sample": adata.obs["sample"].astype(str).values,
        "dataset": adata.obs["dataset"].astype(str).values,
        argmax_col: argmax_int,
    })
    for j, c in enumerate(classes):
        per_cell[f"P_F{c}"] = proba[:, j].astype(np.float32)
    per_cell.to_csv(OUT_PERCELL, sep="\t", index=False,
                    compression="gzip", float_format="%.5f")
    log.info(f"  wrote {OUT_PERCELL} ({len(per_cell):,} rows)")

    # ── Donor-level aggregation + external validation ──────────────────────
    p_cols = [f"P_F{c}" for c in classes]
    donor_agg = aggregate_donor(per_cell, p_cols, argmax_col)
    donor_agg.to_csv(OUT_DONOR_TSV, sep="\t", index=False, float_format="%.6f")
    log.info(f"  wrote {OUT_DONOR_TSV} ({len(donor_agg):,} donors)")

    ext = external_validation(donor_agg, meta, argmax_col)

    # ── Method metrics ─────────────────────────────────────────────────────
    qwk = cv_result["qwk"] if cv_result else float("nan")
    acc = cv_result["accuracy"] if cv_result else float("nan")
    off_by_one = ""
    if cv_result is not None:
        cvp = cv_result["cv"]
        diff = (cvp["F_stage_documented"].astype(int)
                - cvp[f"{argmax_col}_argmax"].astype(int)).abs()
        off_by_one = round(float((diff <= 1).mean()), 4)
    n_train_donors = int(fstage["F_stage_documented"].notna().sum())
    method_row = dict(
        method=f"scanvi_{variant}_v2",
        n_train_donors=n_train_donors,
        n_labelled_cells=n_labelled,
        qwk_cv=round(qwk, 4) if cv_result else "",
        accuracy_cv=round(acc, 4) if cv_result else "",
        off_by_one_cv=off_by_one,
        external_rho=round(ext["rho_overall"], 4),
        n_healthy_F4_misassign=ext["n_healthy_f4"],
        cv_strategy=f"{args.cv_folds}-fold stratified" if cv_result else "skip",
        n_cells=adata.n_obs,
    )
    pd.DataFrame([method_row]).to_csv(OUT_METRICS, sep="\t", index=False)
    log.info(f"  wrote {OUT_METRICS}")

    # ── External validation summary ────────────────────────────────────────
    with open(OUT_SUMMARY, "w") as f:
        f.write(f"# scANVI v2 ({variant}) F-stage external validation\n")
        f.write(f"# cells: {adata.n_obs:,}; n_train donors: {n_train_donors}; "
                f"n_labelled cells: {n_labelled:,}\n")
        f.write(f"# CV: {args.cv_folds}-fold stratified across Andrews documented donors\n")
        if cv_result is not None:
            f.write(f"\nCV QWK    = {cv_result['qwk']:.3f}\n")
            f.write(f"CV accuracy = {cv_result['accuracy']:.3f}\n")
            f.write(f"CV off-by-one tolerance = {off_by_one}\n")
            f.write("\nConfusion matrix (rows=true, cols=pred, F0..F4):\n")
            f.write(str(np.array(cv_result["cm"])) + "\n")
        f.write(f"\nExternal Spearman rho({argmax_col}, disease_stage_numeric)"
                f" = {ext['rho_overall']:.3f} (p={ext['p_overall']:.2e})\n")
        f.write(f"Healthy donors predicted F4: {ext['n_healthy_f4']}\n")
        f.write("\nPer-dataset Spearman rho:\n")
        f.write(ext["per_dataset"].to_string(index=False) + "\n")
        f.write("\nPredicted F-stage x disease_stage_coarse:\n")
        f.write(ext["crosstab"].to_string() + "\n")
    log.info(f"  wrote {OUT_SUMMARY}")

    # ── Training log ───────────────────────────────────────────────────────
    train_record = {
        "script": "03_train_scanvi_v2.py",
        "variant": variant,
        "atlas_in": str(ATLAS_IN),
        "scvi_base": str(SCVI_MODEL_DIR),
        "scanvi_model_dir": str(MODEL_OUT),
        "n_cells": int(adata.n_obs),
        "n_genes_used": int(adata.n_vars),
        "n_labelled_cells": n_labelled,
        "n_train_donors": n_train_donors,
        "hyperparams": {
            "max_epochs": SCANVI_MAX_EPOCHS,
            "n_samples_per_label": SCANVI_N_SAMPLES_PER_LABEL,
            "cv_n_splits": args.cv_folds,
        },
        "cv_qwk": cv_result["qwk"] if cv_result else None,
        "cv_accuracy": cv_result["accuracy"] if cv_result else None,
        "external_rho": ext["rho_overall"],
        "n_healthy_F4_misassign": ext["n_healthy_f4"],
        "scvi_tools_version": scvi.__version__,
        "torch_version": torch.__version__,
        "cuda_available": bool(torch.cuda.is_available()),
        "device_name": (
            torch.cuda.get_device_name(0) if torch.cuda.is_available() else "cpu"
        ),
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    with open(OUT_LOG, "w") as f:
        json.dump(train_record, f, indent=2)
    log.info(f"wrote training log -> {OUT_LOG}")
    log.info("DONE")


if __name__ == "__main__":
    main()
