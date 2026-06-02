#!/usr/bin/env python
"""
343h_scanvi_fstage.py

Train scANVI (semi-supervised scVI) with cell-level F-stage labels broadcast
from donor metadata, then evaluate F-stage prediction accuracy.

Rationale: vanilla scVI's latent is unsupervised — F-stage may not be a
primary axis the model encodes. The 343b/343c/343d pipeline showed donor
pseudobulk + ordinal logistic gives LOOCV QWK 0.742 within Andrews but
external Spearman rho drops to 0.391 and hallucinates F4 in healthy donors
of Liver_Atlas / GSE185477.

scANVI extends scVI with a semi-supervised cell-level label classifier that
forces the latent to capture the supplied labels. By broadcasting each
Andrews donor's F-stage to its hepatocytes (and marking all cells from
non-Andrews donors 'Unknown'), scANVI should learn an F-stage-resolved
latent that transfers better.

Inputs
------
- Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_atlas_annotated.h5ad
    657,804 hepatocytes; .layers['counts'] holds raw counts; HVG-subset already
    applied (2,931 var); _scvi_batch present from the original scVI run.
- Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_documented.tsv
    58 Andrews donors with documented Kleiner F-stage (F0=8, F1=8, F2=12, F3=11, F4=19).
- Persisted scVI model at:
    Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/scvi_hepatocyte_model/model.pt

Outputs
-------
- donor_fstage_scanvi_predicted.tsv
    sample, dataset, F_stage_scanvi_argmax, F_stage_scanvi_mode,
    posterior_F0..F4_mean, n_hepatocytes
- percell_fstage_scanvi.tsv.gz
    cell_id, sample, dataset, F_stage_scanvi, P_F0..P_F4
- scanvi_method_metrics.tsv
- Append row to fstage_method_comparison.tsv
- Models: results_gpu_v2/ccc/stage_trajectory/scanvi_model/{scvi_full,scanvi_full}

CV strategy
-----------
5-fold stratified CV across the 58 Andrews donors (instead of full LOOCV)
because retraining scANVI 58 times is prohibitive. Each fold trains a new
scANVI model on the four non-held-out folds (with held-out donors' cells
relabelled 'Unknown'), then predicts the held-out donors' F-stage. Discloses
this trade-off in the metrics tsv.

Implementation
--------------
- Loads the persisted scVI model and reattaches the AnnData (HVG-subset).
- Subsets to hepatocytes (atlas is already hepatocytes only).
- Optional dataset-stratified subsample if GPU memory pressure arises.
- scVI is already trained; we DO NOT retrain it. scANVI is initialised from
  this scVI via SCANVI.from_scvi_model(...) and only the classifier head + a
  short fine-tune are trained.
"""
from __future__ import annotations
import os
import sys
import time
import json
import gc
import argparse
from pathlib import Path
from collections import Counter

# Disable scvi-tools internet calls and wandb logging
os.environ.setdefault("WANDB_MODE", "disabled")
os.environ.setdefault("WANDB_DISABLED", "true")
os.environ.setdefault("SCVI_NO_TELEMETRY", "true")

import numpy as np
import pandas as pd
import anndata as ad

import torch
import scvi

from sklearn.metrics import cohen_kappa_score, confusion_matrix, accuracy_score
from sklearn.model_selection import StratifiedKFold
from scipy import stats


# ── Paths ────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))

HEP_H5AD = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_atlas_annotated.h5ad"
SCVI_MODEL_DIR = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/scvi_hepatocyte_model"
FSTAGE_DOC = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_documented.tsv"
META_EXT = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"

OUT_DIR = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
MODEL_OUT = OUT_DIR / "scanvi_model"
OUT_DONOR_TSV = OUT_DIR / "donor_fstage_scanvi_predicted.tsv"
OUT_PERCELL = OUT_DIR / "percell_fstage_scanvi.tsv.gz"
OUT_METRICS = OUT_DIR / "scanvi_method_metrics.tsv"
OUT_COMPARE = OUT_DIR / "fstage_method_comparison.tsv"

# scVI / scANVI training hyperparameters (match original scVI run)
SCVI_N_LATENT = 20
SCVI_N_LAYERS = 2
SCVI_N_HIDDEN = 128
SCVI_GENE_LIKELIHOOD = "nb"
SCVI_DROPOUT = 0.1

# scANVI fine-tune
SCANVI_MAX_EPOCHS = 75
SCANVI_N_SAMPLES_PER_LABEL = 100

UNLABELED = "Unknown"
F_STAGES = [0, 1, 2, 3, 4]
F_STAGE_STR = [str(f) for f in F_STAGES]


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def set_seeds(seed: int = 42):
    import random
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed); np.random.seed(seed)
    torch.manual_seed(seed); torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    scvi.settings.seed = seed


def broadcast_labels(adata: ad.AnnData, fstage_df: pd.DataFrame) -> pd.Series:
    """Cell-level fstage label series; documented donors get '0'..'4', else 'Unknown'."""
    fstage_df = fstage_df.dropna(subset=["F_stage_documented"]).copy()
    fstage_df["F_stage_documented"] = pd.to_numeric(fstage_df["F_stage_documented"], errors="coerce")
    fstage_df = fstage_df.dropna(subset=["F_stage_documented"])
    fstage_df["F_stage_documented"] = fstage_df["F_stage_documented"].astype(int)

    donor_to_f = dict(zip(fstage_df["sample"], fstage_df["F_stage_documented"]))
    labels = adata.obs["sample"].map(donor_to_f).astype("object")
    labels = labels.where(labels.notna(), UNLABELED)
    labels = labels.map(lambda v: str(int(v)) if v != UNLABELED else UNLABELED)
    n_lab = (labels != UNLABELED).sum()
    log(f"  broadcast: {n_lab:,}/{len(labels):,} cells labelled (Andrews donors with F-stage)")
    log(f"  per-class counts: {dict(Counter(labels))}")
    return labels


def _register_null_io_fallback():
    """The hepatocyte atlas was written with anndata>=0.10 using encoding_type='null'
    for uns/log1p/base. The reader version in scanpy_env doesn't yet have a
    registered handler. Register a no-op (returns None) fallback so the file loads.
    """
    try:
        import h5py
        from anndata._io.specs.registry import _REGISTRY, IOSpec
        _REGISTRY.register_read(h5py.Dataset, IOSpec("null", "0.1.0"))(lambda elem, _reader=None: None)
        _REGISTRY.register_read(h5py.Group,  IOSpec("null", "0.1.0"))(lambda elem, _reader=None: None)
    except Exception as e:
        log(f"  WARN: failed to register null IO fallback ({e})")


def load_atlas_with_counts() -> ad.AnnData:
    _register_null_io_fallback()
    log(f"loading {HEP_H5AD.name}")
    a = ad.read_h5ad(HEP_H5AD)
    log(f"  shape: {a.shape}; layers: {list(a.layers.keys())}; obsm: {list(a.obsm.keys())}")
    if "counts" not in a.layers:
        sys.exit("FATAL: 'counts' layer missing")
    # Use counts as X for scvi (setup_anndata layer arg)
    return a


def make_scanvi(adata: ad.AnnData, labels_key: str, scvi_model=None) -> "scvi.model.SCANVI":
    """Construct scANVI: load pre-trained scVI then convert."""
    # If a scvi_model is supplied, transfer; otherwise train scVI from scratch.
    if scvi_model is None:
        log("  training scVI from scratch (no persisted model found / chose to retrain)")
        scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="dataset")
        scvi_model = scvi.model.SCVI(
            adata,
            n_latent=SCVI_N_LATENT, n_layers=SCVI_N_LAYERS, n_hidden=SCVI_N_HIDDEN,
            gene_likelihood=SCVI_GENE_LIKELIHOOD, dropout_rate=SCVI_DROPOUT,
        )
        scvi_model.train(
            max_epochs=150, early_stopping=True, early_stopping_patience=10,
            train_size=0.9, batch_size=512, enable_progress_bar=False,
        )
    # Note: scANVI needs the SAME AnnData object as scVI was trained on (or
    # a subset). For LOOCV we mutate adata.obs[labels_key] but keep cell order
    # stable.
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model, labels_key=labels_key, unlabeled_category=UNLABELED,
        adata=adata,
    )
    scanvi_model.train(
        max_epochs=SCANVI_MAX_EPOCHS,
        n_samples_per_label=SCANVI_N_SAMPLES_PER_LABEL,
        check_val_every_n_epoch=10,
        enable_progress_bar=False,
    )
    return scanvi_model


def aggregate_donor(per_cell: pd.DataFrame, p_cols: list) -> pd.DataFrame:
    """Donor-level aggregation: argmax of mean posterior + mode of per-cell argmax."""
    grouped_mean = per_cell.groupby(["sample", "dataset"], observed=True)[p_cols].mean()
    grouped_mean.columns = [f"posterior_F{c.split('P_F')[-1]}_mean" for c in p_cols]
    argmax_per_donor = grouped_mean.values.argmax(axis=1)
    mode_per_donor = (per_cell.groupby(["sample", "dataset"], observed=True)["F_stage_scanvi"]
                              .agg(lambda x: x.value_counts().idxmax()))
    n_hep = per_cell.groupby(["sample", "dataset"], observed=True).size().rename("n_hepatocytes")
    out = grouped_mean.copy()
    out["F_stage_scanvi_argmax"] = argmax_per_donor
    out["F_stage_scanvi_mode"] = mode_per_donor.values
    out["n_hepatocytes"] = n_hep.values
    return out.reset_index()


def predict_with_proba(scanvi_model, adata) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Return (argmax_int, P[N,5_in_F-order], class_labels_str)."""
    soft = scanvi_model.predict(adata, soft=True)  # DataFrame: index=cells, cols=labels
    # remove Unknown column if present
    if UNLABELED in soft.columns:
        soft = soft.drop(columns=[UNLABELED])
    # reindex to consistent F0..F4 order
    for s in F_STAGE_STR:
        if s not in soft.columns:
            soft[s] = 0.0
    soft = soft[F_STAGE_STR]
    rowsum = soft.values.sum(axis=1, keepdims=True)
    rowsum[rowsum == 0] = 1.0
    proba = soft.values / rowsum
    argmax = proba.argmax(axis=1)
    return argmax, proba, F_STAGE_STR


def fold_cv(adata: ad.AnnData, donor_labels: pd.DataFrame, scvi_model_main,
            n_splits: int = 5, seed: int = 42) -> dict:
    """5-fold stratified CV across Andrews donors. Returns dict of QWK / accuracy / per-donor preds."""
    docs = donor_labels.dropna(subset=["F_stage_documented"]).copy()
    docs["F_stage_documented"] = docs["F_stage_documented"].astype(int)
    donors = docs["sample"].values
    y = docs["F_stage_documented"].values

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    fold_preds = []
    for fold_idx, (train_idx, test_idx) in enumerate(skf.split(donors, y)):
        held_out_donors = set(donors[test_idx])
        log(f"FOLD {fold_idx+1}/{n_splits}: holding out {len(held_out_donors)} donors")
        # Build a fold-specific labels column: held-out donors marked UNLABELED
        cell_donor = adata.obs["sample"].astype(str).values
        full_labels = broadcast_labels(adata, docs).values
        fold_labels = full_labels.copy()
        for i, s in enumerate(cell_donor):
            if s in held_out_donors:
                fold_labels[i] = UNLABELED
        col = f"fstage_fold{fold_idx}"
        adata.obs[col] = pd.Categorical(fold_labels, categories=F_STAGE_STR + [UNLABELED])

        # SCANVI needs the labels_key to be a column on the same AnnData used by SCVI
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
        # Predict on held-out cells
        held_mask = adata.obs["sample"].isin(held_out_donors).values
        a_held = adata[held_mask]
        argmax_int, proba, classes = predict_with_proba(scanvi, a_held)
        per_cell = pd.DataFrame({
            "sample":  a_held.obs["sample"].values,
            "dataset": a_held.obs["dataset"].values,
            "F_stage_scanvi": argmax_int,
        })
        for j, c in enumerate(classes):
            per_cell[f"P_F{c}"] = proba[:, j]
        p_cols = [f"P_F{c}" for c in classes]
        donor_agg = aggregate_donor(per_cell, p_cols)
        donor_agg["fold"] = fold_idx
        fold_preds.append(donor_agg)

        # cleanup
        del scanvi
        gc.collect()
        torch.cuda.empty_cache() if torch.cuda.is_available() else None

    cv = pd.concat(fold_preds, ignore_index=True)
    cv = cv.merge(docs[["sample", "F_stage_documented"]], on="sample", how="left")
    qwk = cohen_kappa_score(cv["F_stage_documented"].astype(int),
                            cv["F_stage_scanvi_argmax"].astype(int),
                            weights="quadratic", labels=F_STAGES)
    acc = accuracy_score(cv["F_stage_documented"].astype(int),
                         cv["F_stage_scanvi_argmax"].astype(int))
    cm = confusion_matrix(cv["F_stage_documented"].astype(int),
                          cv["F_stage_scanvi_argmax"].astype(int),
                          labels=F_STAGES)
    log(f"5-fold CV: QWK = {qwk:.3f}, accuracy = {acc:.3f}")
    log("Confusion matrix (rows=true, cols=pred):")
    log(str(cm))
    return dict(qwk=float(qwk), accuracy=float(acc), cm=cm.tolist(), cv=cv, n_splits=n_splits)


def external_validation(donor_preds: pd.DataFrame, meta: pd.DataFrame) -> dict:
    """Spearman rho(predicted, disease_stage_numeric) overall and per-dataset; Healthy->F4 errors."""
    df = donor_preds.merge(meta[["sample", "disease_stage_coarse", "disease_stage_numeric"]],
                            on="sample", how="left")
    df = df.dropna(subset=["F_stage_scanvi_argmax", "disease_stage_numeric"])
    rho_overall, p_overall = stats.spearmanr(df["F_stage_scanvi_argmax"], df["disease_stage_numeric"])
    log(f"External Spearman rho(F_stage_scanvi, disease_stage_numeric) = {rho_overall:.3f} (n={len(df)})")
    per_ds = []
    for ds, g in df.groupby("dataset"):
        if g["disease_stage_numeric"].nunique() < 2:
            per_ds.append(dict(dataset=ds, rho=np.nan, p=np.nan, n=len(g)))
            continue
        r, p = stats.spearmanr(g["F_stage_scanvi_argmax"], g["disease_stage_numeric"])
        per_ds.append(dict(dataset=ds, rho=float(r), p=float(p), n=int(len(g))))
    per_ds = pd.DataFrame(per_ds)
    n_healthy_f4 = int(((df["F_stage_scanvi_argmax"] == 4) &
                        (df["disease_stage_coarse"] == "Healthy")).sum())
    cross = pd.crosstab(df["F_stage_scanvi_argmax"], df["disease_stage_coarse"]).reindex(
        index=F_STAGES, columns=["Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"],
        fill_value=0,
    )
    log(f"Healthy donors predicted as F4: {n_healthy_f4}")
    log("Predicted F-stage x disease_stage_coarse:")
    log(str(cross))
    return dict(
        rho_overall=float(rho_overall), p_overall=float(p_overall),
        per_dataset=per_ds, n_healthy_f4=n_healthy_f4, crosstab=cross,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subsample", type=int, default=0,
                    help="If >0, subsample atlas to this many cells (stratified by dataset).")
    ap.add_argument("--retrain-scvi", action="store_true",
                    help="If set, retrain scVI from scratch (default: load persisted scVI).")
    ap.add_argument("--skip-cv", action="store_true",
                    help="Skip 5-fold CV (faster: just train once and predict).")
    ap.add_argument("--cv-folds", type=int, default=5)
    args = ap.parse_args()

    set_seeds(42)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    MODEL_OUT.mkdir(parents=True, exist_ok=True)

    log(f"GPU available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        log(f"  device: {torch.cuda.get_device_name(0)}")

    adata = load_atlas_with_counts()

    if args.subsample and args.subsample < adata.n_obs:
        rng = np.random.default_rng(42)
        idx = (
            pd.Series(range(adata.n_obs))
            .groupby(adata.obs["dataset"].values, observed=True)
            .sample(frac=args.subsample / adata.n_obs, random_state=42)
            .index.values
        )
        idx = np.array(sorted(set(idx)))
        adata = adata[idx].copy()
        log(f"subsampled to {adata.n_obs:,} cells (stratified by dataset)")

    fstage = pd.read_csv(FSTAGE_DOC, sep="\t")
    fstage["F_stage_documented"] = pd.to_numeric(fstage["F_stage_documented"], errors="coerce")
    meta = pd.read_csv(META_EXT, sep="\t")

    # Cell-level labels
    cell_labels = broadcast_labels(adata, fstage)
    adata.obs["fstage_label"] = pd.Categorical(cell_labels, categories=F_STAGE_STR + [UNLABELED])

    # ── Load (or train) base scVI ─────────────────────────────────────────
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="dataset")
    if (not args.retrain_scvi) and (SCVI_MODEL_DIR / "model.pt").exists():
        try:
            log(f"loading persisted scVI from {SCVI_MODEL_DIR}")
            scvi_model = scvi.model.SCVI.load(str(SCVI_MODEL_DIR), adata=adata)
            log("  loaded")
        except Exception as e:
            log(f"  WARNING: failed to load persisted scVI ({e}); training fresh")
            scvi_model = scvi.model.SCVI(
                adata,
                n_latent=SCVI_N_LATENT, n_layers=SCVI_N_LAYERS, n_hidden=SCVI_N_HIDDEN,
                gene_likelihood=SCVI_GENE_LIKELIHOOD, dropout_rate=SCVI_DROPOUT,
            )
            scvi_model.train(
                max_epochs=150, early_stopping=True, early_stopping_patience=10,
                train_size=0.9, batch_size=512, enable_progress_bar=False,
            )
            scvi_model.save(str(MODEL_OUT / "scvi_full"), overwrite=True)
    else:
        log("training scVI from scratch (retrain-scvi flag or model missing)")
        scvi_model = scvi.model.SCVI(
            adata,
            n_latent=SCVI_N_LATENT, n_layers=SCVI_N_LAYERS, n_hidden=SCVI_N_HIDDEN,
            gene_likelihood=SCVI_GENE_LIKELIHOOD, dropout_rate=SCVI_DROPOUT,
        )
        scvi_model.train(
            max_epochs=150, early_stopping=True, early_stopping_patience=10,
            train_size=0.9, batch_size=512, enable_progress_bar=False,
        )
        scvi_model.save(str(MODEL_OUT / "scvi_full"), overwrite=True)

    # ── 5-fold CV (within Andrews) ────────────────────────────────────────
    cv_result = None
    if not args.skip_cv:
        log("===== 5-fold CV =====")
        try:
            cv_result = fold_cv(adata, fstage, scvi_model, n_splits=args.cv_folds, seed=42)
        except Exception as e:
            log(f"CV FAILED: {e}")
            import traceback; traceback.print_exc()
            cv_result = None

    # ── Final scANVI: train on full Andrews labels and predict everything ──
    log("===== final scANVI on full Andrews labels =====")
    # Reset labels to the full broadcast
    adata.obs["fstage_label"] = pd.Categorical(cell_labels, categories=F_STAGE_STR + [UNLABELED])
    scanvi_full = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        labels_key="fstage_label",
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
    log(f"  saved scANVI model -> {MODEL_OUT / 'scanvi_full'}")

    # ── Predict for every hepatocyte ──────────────────────────────────────
    log("predicting F-stage for all hepatocytes")
    argmax_int, proba, classes = predict_with_proba(scanvi_full, adata)
    per_cell = pd.DataFrame({
        "cell_id": adata.obs_names.values.astype(str),
        "sample":  adata.obs["sample"].astype(str).values,
        "dataset": adata.obs["dataset"].astype(str).values,
        "F_stage_scanvi": argmax_int,
    })
    for j, c in enumerate(classes):
        per_cell[f"P_F{c}"] = proba[:, j].astype(np.float32)

    # write per-cell (compressed)
    per_cell.to_csv(OUT_PERCELL, sep="\t", index=False, compression="gzip", float_format="%.5f")
    log(f"  wrote {OUT_PERCELL} ({len(per_cell):,} rows)")

    # ── Donor-level aggregation ────────────────────────────────────────────
    p_cols = [f"P_F{c}" for c in classes]
    donor_agg = aggregate_donor(per_cell, p_cols)
    donor_agg.to_csv(OUT_DONOR_TSV, sep="\t", index=False, float_format="%.6f")
    log(f"  wrote {OUT_DONOR_TSV} ({len(donor_agg):,} donors)")

    # ── External validation ──────────────────────────────────────────────
    ext = external_validation(donor_agg, meta)

    # ── Write method metrics + append to comparison ──────────────────────
    qwk = cv_result["qwk"] if cv_result else float("nan")
    acc = cv_result["accuracy"] if cv_result else float("nan")
    n_train = int(fstage["F_stage_documented"].notna().sum())
    method_row = dict(
        method="scanvi",
        n_train=n_train,
        qwk_loocv=round(qwk, 4) if cv_result else "",
        accuracy=round(acc, 4) if cv_result else "",
        external_rho=round(ext["rho_overall"], 4),
        n_healthy_F4_misassign=ext["n_healthy_f4"],
        cv_strategy=f"{args.cv_folds}-fold stratified" if cv_result else "skip",
        n_cells=adata.n_obs,
    )
    pd.DataFrame([method_row]).to_csv(OUT_METRICS, sep="\t", index=False)
    log(f"  wrote {OUT_METRICS}")

    # Append to fstage_method_comparison.tsv (preserve existing rows)
    if OUT_COMPARE.exists():
        comp = pd.read_csv(OUT_COMPARE, sep="\t")
    else:
        comp = pd.DataFrame(columns=["method", "n_train", "qwk_loocv", "accuracy_loocv", "off_by_one_loocv"])
    # remove any prior scanvi row to avoid duplicates
    comp = comp[comp["method"] != "scanvi"].copy()
    # compute off-by-one tolerance from CV preds
    off_by_one = ""
    if cv_result is not None:
        cv_preds = cv_result["cv"]
        diff = (cv_preds["F_stage_documented"].astype(int) - cv_preds["F_stage_scanvi_argmax"].astype(int)).abs()
        off_by_one = round(float((diff <= 1).mean()), 4)
    new_row = dict(
        method="scanvi",
        n_train=n_train,
        qwk_loocv=round(qwk, 4) if cv_result else "",
        accuracy_loocv=round(acc, 4) if cv_result else "",
        off_by_one_loocv=off_by_one,
    )
    comp = pd.concat([comp, pd.DataFrame([new_row])], ignore_index=True)
    comp.to_csv(OUT_COMPARE, sep="\t", index=False)
    log(f"  appended scanvi row -> {OUT_COMPARE}")

    # ── Write external-validation summary ────────────────────────────────
    summary_path = OUT_DIR / "scanvi_external_validation_summary.txt"
    with open(summary_path, "w") as f:
        f.write("# scANVI F-stage external validation\n")
        f.write(f"# cells: {adata.n_obs:,}; n_train donors: {n_train}\n")
        f.write(f"# CV: {args.cv_folds}-fold stratified across Andrews donors (not LOOCV because retraining 58x is prohibitive)\n")
        if cv_result is not None:
            f.write(f"\nCV QWK    = {cv_result['qwk']:.3f}\n")
            f.write(f"CV accuracy = {cv_result['accuracy']:.3f}\n")
            f.write(f"CV off-by-one tolerance = {off_by_one}\n")
            f.write("\nConfusion matrix (rows=true, cols=pred, F0..F4):\n")
            f.write(str(np.array(cv_result["cm"])) + "\n")
        f.write(f"\nExternal Spearman rho(F_stage_scanvi, disease_stage_numeric) = {ext['rho_overall']:.3f} (p={ext['p_overall']:.2e})\n")
        f.write(f"Healthy donors predicted as F4: {ext['n_healthy_f4']}\n")
        f.write("\nPer-dataset Spearman rho:\n")
        f.write(ext["per_dataset"].to_string(index=False) + "\n")
        f.write("\nPredicted F-stage x disease_stage_coarse:\n")
        f.write(ext["crosstab"].to_string() + "\n")
        # head-to-head verdict
        BASELINE_RHO = 0.391
        BASELINE_F4ERR = 27
        better_rho = ext["rho_overall"] > BASELINE_RHO
        better_err = ext["n_healthy_f4"] < BASELINE_F4ERR
        f.write("\nHead-to-head vs vanilla scVI ordinal:\n")
        f.write(f"  external rho:    scANVI={ext['rho_overall']:.3f} vs vanilla=0.391  -> {'BETTER' if better_rho else 'NOT BETTER'}\n")
        f.write(f"  Healthy->F4:     scANVI={ext['n_healthy_f4']}  vs vanilla=27       -> {'BETTER' if better_err else 'NOT BETTER'}\n")
    log(f"  wrote {summary_path}")

    log("DONE")


if __name__ == "__main__":
    main()
