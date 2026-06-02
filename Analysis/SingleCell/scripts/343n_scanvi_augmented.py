#!/usr/bin/env python
"""
343n_scanvi_augmented.py

Train scANVI with AUGMENTED cell-level F-stage labels:
  - Andrews documented donor F-stage      -> cell label = F-stage (0..4)
  - disease_stage_coarse == "Healthy" AND not dubious_healthy -> cell label = 0 (F0)
  - disease_stage_coarse == "Cirrhosis"                       -> cell label = 4 (F4)
  - otherwise                                                  -> cell label = "Unknown"

This extends the cell-level scANVI framework (Script 343h, vanilla Andrews-only)
with the donor-level augmentation trick that gave Script 343m the external
Spearman rho jump from 0.391 (vanilla scVI ordinal) to 0.639 (augmented scVI
ordinal).

Strategy
--------
- Load persisted scVI base model (frozen).
- Build cell-level labels with the broadcast above (~515K labelled).
- Train scANVI head with n_samples_per_label=200 (to balance large F0/F4
  anchor sets against smaller F1-F3 documented sets).
- 5-fold CV across Andrews 58 documented donors. In each fold, hold out the
  fold's Andrews donors' cells (mark "Unknown") while keeping ALL clean-healthy
  and cirrhosis anchors in training, then refit scANVI head from the frozen
  scVI. (Matches the 343h CV strategy; full LOOCV across 58 documented donors
  is too expensive given each scANVI head fit is ~1h on L40S.)
- Predict F-stage per cell via scanvi_model.predict(soft=True).
- Aggregate to per-donor argmax + posterior-mean.
- Evaluate external Spearman rho vs disease_stage_numeric and Healthy->F4 errors.

Outputs (Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/)
-----------------------------------------------------------------
- donor_fstage_scanvi_augmented_predicted.tsv
- percell_fstage_scanvi_augmented.tsv.gz
- scanvi_augmented_method_metrics.tsv
- appended row in fstage_method_comparison.tsv
- scanvi_augmented_external_validation_summary.txt
- figures/supplementary/stage_ccc/figS_scanvi_augmented_comparison.pdf
- Models: scanvi_augmented_model/scanvi_full/
"""
from __future__ import annotations
import os
import sys
import time
import gc
import argparse
from pathlib import Path
from collections import Counter

# Disable scvi-tools internet calls + wandb logging
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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# ── Paths ────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))

HEP_H5AD = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_atlas_annotated.h5ad"
SCVI_MODEL_DIR = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/scvi_hepatocyte_model"
STAGE_DIR = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
FSTAGE_DOC = STAGE_DIR / "donor_fstage_documented.tsv"
DUBIOUS = STAGE_DIR / "dubious_healthy_donors.tsv"
META_EXT = STAGE_DIR / "donor_metadata_extended.tsv"

OUT_DIR = STAGE_DIR
MODEL_OUT = OUT_DIR / "scanvi_augmented_model"
OUT_DONOR_TSV = OUT_DIR / "donor_fstage_scanvi_augmented_predicted.tsv"
OUT_PERCELL = OUT_DIR / "percell_fstage_scanvi_augmented.tsv.gz"
OUT_METRICS = OUT_DIR / "scanvi_augmented_method_metrics.tsv"
OUT_COMPARE = OUT_DIR / "fstage_method_comparison.tsv"
OUT_SUMMARY = OUT_DIR / "scanvi_augmented_external_validation_summary.txt"
OUT_PDF = PROJECT_ROOT / "figures/supplementary/stage_ccc/figS_scanvi_augmented_comparison.pdf"

# scVI / scANVI training hyperparameters (match Script 343h)
SCVI_N_LATENT = 20
SCVI_N_LAYERS = 2
SCVI_N_HIDDEN = 128
SCVI_GENE_LIKELIHOOD = "nb"
SCVI_DROPOUT = 0.1

# scANVI fine-tune: n_samples_per_label=200 (vs 100 in 343h vanilla) to balance
# large F0/F4 anchor sets against smaller F1-F3 documented sets.
SCANVI_MAX_EPOCHS = 75
SCANVI_N_SAMPLES_PER_LABEL = 200

UNLABELED = "Unknown"
F_STAGES = [0, 1, 2, 3, 4]
F_STAGE_STR = [str(f) for f in F_STAGES]


MASLD_PAL = {
    "Healthy":         "#9E9E9E",
    "Steatosis":       "#FFB454",
    "Steatohepatitis": "#E66F51",
    "Cirrhosis":       "#7A0E0E",
    "F0":              "#9E9E9E",
    "F1":              "#FFD16B",
    "F2":              "#FFB454",
    "F3":              "#E66F51",
    "F4":              "#7A0E0E",
    "vanilla_ord":     "#5B7FB1",
    "scanvi_vanilla":  "#8A6FB3",
    "augmented_ord":   "#2E8B57",
    "scanvi_aug":      "#1B6B3A",
}


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


def _register_null_io_fallback():
    """anndata>=0.10 may write encoding_type='null' for uns/log1p/base; older
    readers crash on this. Register a no-op."""
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
    return a


def broadcast_augmented_labels(
    adata: ad.AnnData,
    fstage_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    dubious_samples: set,
) -> pd.Series:
    """Cell-level fstage label series with augmentation.

    Per-cell rule (in priority order):
      1. donor in Andrews documented      -> F-stage 0..4
      2. donor coarse == 'Healthy' and not dubious -> 0 (F0)
      3. donor coarse == 'Cirrhosis'              -> 4 (F4)
      4. else                                      -> 'Unknown'
    """
    # documented donors
    fstage_df = fstage_df.dropna(subset=["F_stage_documented"]).copy()
    fstage_df["F_stage_documented"] = pd.to_numeric(
        fstage_df["F_stage_documented"], errors="coerce")
    fstage_df = fstage_df.dropna(subset=["F_stage_documented"])
    fstage_df["F_stage_documented"] = fstage_df["F_stage_documented"].astype(int)
    docs = dict(zip(fstage_df["sample"], fstage_df["F_stage_documented"]))

    # coarse stage
    coarse = dict(zip(meta_df["sample"], meta_df["disease_stage_coarse"]))

    samples = adata.obs["sample"].astype(str).values
    labels = np.full(len(samples), UNLABELED, dtype=object)
    n_doc = 0
    n_chealthy = 0
    n_cirr = 0
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
    n_lab = (labels != UNLABELED).sum()
    log(f"  broadcast augmented: {n_lab:,}/{len(labels):,} cells labelled")
    log(f"    Andrews documented   : {n_doc:,}")
    log(f"    clean-healthy -> F0  : {n_chealthy:,}")
    log(f"    Cirrhosis -> F4      : {n_cirr:,}")
    log(f"    per-class counts: {dict(Counter(labels))}")
    return labels


def predict_with_proba(scanvi_model, adata) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Return (argmax_int, P[N,5_in_F-order], class_labels_str)."""
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


def aggregate_donor(per_cell: pd.DataFrame, p_cols: list) -> pd.DataFrame:
    """Donor-level aggregation: argmax of mean posterior + mode of per-cell argmax."""
    grouped_mean = per_cell.groupby(["sample", "dataset"], observed=True)[p_cols].mean()
    grouped_mean.columns = [f"posterior_F{c.split('P_F')[-1]}_mean" for c in p_cols]
    argmax_per_donor = grouped_mean.values.argmax(axis=1)
    mode_per_donor = (per_cell.groupby(["sample", "dataset"], observed=True)["F_stage_scanvi_aug"]
                              .agg(lambda x: x.value_counts().idxmax()))
    n_hep = per_cell.groupby(["sample", "dataset"], observed=True).size().rename("n_hepatocytes")
    posterior_mean = (
        grouped_mean.values * np.array([0, 1, 2, 3, 4])
    ).sum(axis=1)
    out = grouped_mean.copy()
    out["F_stage_scanvi_augmented_argmax"] = argmax_per_donor
    out["F_stage_scanvi_augmented_mode"] = mode_per_donor.values
    out["posterior_mean"] = posterior_mean
    out["n_cells"] = n_hep.values
    return out.reset_index()


def fold_cv(
    adata: ad.AnnData,
    fstage_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    dubious_samples: set,
    full_cell_labels: pd.Series,
    scvi_model_main,
    n_splits: int = 5,
    seed: int = 42,
) -> dict:
    """K-fold stratified CV across Andrews documented donors.

    For each fold: held-out Andrews donors' cells are relabelled 'Unknown'.
    All clean-healthy and cirrhosis anchors stay in training. Refit scANVI
    head from the frozen scVI base, predict the held-out cells, aggregate to
    donor level (argmax of mean posterior).
    """
    docs = fstage_df.dropna(subset=["F_stage_documented"]).copy()
    docs["F_stage_documented"] = docs["F_stage_documented"].astype(int)
    donors = docs["sample"].values
    y = docs["F_stage_documented"].values

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    fold_preds = []
    full_cell_labels_arr = full_cell_labels.values.copy()
    cell_samples = adata.obs["sample"].astype(str).values

    for fold_idx, (train_idx, test_idx) in enumerate(skf.split(donors, y)):
        held_out_donors = set(donors[test_idx])
        log(f"FOLD {fold_idx+1}/{n_splits}: holding out {len(held_out_donors)} "
            f"Andrews donors (clean-healthy & cirrhosis anchors stay in training)")

        # Build fold-specific labels: held-out Andrews donors' cells -> Unknown
        fold_labels = full_cell_labels_arr.copy()
        held_mask = np.isin(cell_samples, list(held_out_donors))
        fold_labels[held_mask] = UNLABELED

        n_held_cells = int(held_mask.sum())
        n_train_lab = int((fold_labels != UNLABELED).sum())
        log(f"  cells held out: {n_held_cells:,}; "
            f"training-labelled: {n_train_lab:,}; "
            f"per-class counts: {dict(Counter(fold_labels))}")

        col = f"fstage_aug_fold{fold_idx}"
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
            "sample":  a_held.obs["sample"].astype(str).values,
            "dataset": a_held.obs["dataset"].astype(str).values,
            "F_stage_scanvi_aug": argmax_int,
        })
        for j, c in enumerate(classes):
            per_cell[f"P_F{c}"] = proba[:, j]
        p_cols = [f"P_F{c}" for c in classes]
        donor_agg = aggregate_donor(per_cell, p_cols)
        donor_agg["fold"] = fold_idx
        fold_preds.append(donor_agg)

        del scanvi
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    cv = pd.concat(fold_preds, ignore_index=True)
    cv = cv.merge(docs[["sample", "F_stage_documented"]],
                   on="sample", how="left")
    qwk = cohen_kappa_score(
        cv["F_stage_documented"].astype(int),
        cv["F_stage_scanvi_augmented_argmax"].astype(int),
        weights="quadratic", labels=F_STAGES)
    acc = accuracy_score(
        cv["F_stage_documented"].astype(int),
        cv["F_stage_scanvi_augmented_argmax"].astype(int))
    cm = confusion_matrix(
        cv["F_stage_documented"].astype(int),
        cv["F_stage_scanvi_augmented_argmax"].astype(int),
        labels=F_STAGES)
    log(f"{n_splits}-fold CV (Andrews donors): QWK={qwk:.3f}, accuracy={acc:.3f}")
    log("Confusion matrix (rows=true F-stage, cols=pred F-stage):")
    log(str(cm))
    return dict(qwk=float(qwk), accuracy=float(acc), cm=cm.tolist(),
                cv=cv, n_splits=n_splits)


def external_validation(donor_preds: pd.DataFrame, meta: pd.DataFrame) -> dict:
    df = donor_preds.merge(
        meta[["sample", "disease_stage_coarse", "disease_stage_numeric"]],
        on="sample", how="left")
    df = df.dropna(subset=["F_stage_scanvi_augmented_argmax",
                            "disease_stage_numeric"])
    rho_overall, p_overall = stats.spearmanr(
        df["F_stage_scanvi_augmented_argmax"],
        df["disease_stage_numeric"])
    log(f"External Spearman rho(F_stage_scanvi_aug, disease_stage_numeric) "
        f"= {rho_overall:.3f} (n={len(df)})")
    per_ds = []
    for ds, g in df.groupby("dataset"):
        if g["disease_stage_numeric"].nunique() < 2:
            per_ds.append(dict(dataset=ds, rho=np.nan, p=np.nan, n=len(g)))
            continue
        r, p = stats.spearmanr(g["F_stage_scanvi_augmented_argmax"],
                                g["disease_stage_numeric"])
        per_ds.append(dict(dataset=ds, rho=float(r), p=float(p), n=int(len(g))))
    per_ds = pd.DataFrame(per_ds)
    n_healthy_f4 = int(((df["F_stage_scanvi_augmented_argmax"] == 4) &
                        (df["disease_stage_coarse"] == "Healthy")).sum())
    cross = pd.crosstab(
        df["F_stage_scanvi_augmented_argmax"],
        df["disease_stage_coarse"]).reindex(
        index=F_STAGES,
        columns=["Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"],
        fill_value=0,
    )

    # SH breakdown
    sh = df[df["disease_stage_coarse"] == "Steatohepatitis"]
    sh_counts = sh["F_stage_scanvi_augmented_argmax"].value_counts().reindex(
        F_STAGES, fill_value=0)
    sh_frac = (sh_counts / max(len(sh), 1)).round(4)

    log(f"Healthy donors predicted as F4: {n_healthy_f4}")
    log("Predicted F-stage x disease_stage_coarse:")
    log(str(cross))
    log(f"Steatohepatitis breakdown (n={len(sh)}):")
    for k in F_STAGES:
        log(f"  F{k}: {int(sh_counts[k]):3d}  ({sh_frac[k]*100:5.1f}%)")
    return dict(
        rho_overall=float(rho_overall), p_overall=float(p_overall),
        per_dataset=per_ds, n_healthy_f4=n_healthy_f4, crosstab=cross,
        n_sh=int(len(sh)), sh_counts=sh_counts.tolist(),
        sh_frac=sh_frac.tolist(),
    )


def make_figure(cv_result, ext, donor_agg, meta, n_train_cells):
    """4-panel head-to-head with augmented scVI ordinal and vanilla baselines."""
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)

    # Reference numbers (from Script 343m and 343b)
    REF = {
        "vanilla_scvi_ord": dict(qwk=0.7422, rho=0.391, h2f4=27,
                                  sh_frac=[0.069, 0.140, 0.187, 0.314, 0.291]),
        "augmented_scvi_ord": dict(qwk=0.720, rho=0.639, h2f4=2,
                                    sh_frac=[0.00, 0.11, 0.15, 0.43, 0.32]),
    }

    with PdfPages(OUT_PDF) as pdf:
        fig, axes = plt.subplots(2, 2, figsize=(13, 11))

        # Panel A: head-to-head QWK and external rho
        ax = axes[0, 0]
        methods = ["vanilla scVI\nordinal (343b)",
                   "augmented scVI\nordinal (343m)",
                   "scANVI\naugmented (this)"]
        qwks = [REF["vanilla_scvi_ord"]["qwk"],
                REF["augmented_scvi_ord"]["qwk"],
                cv_result["qwk"] if cv_result else np.nan]
        rhos = [REF["vanilla_scvi_ord"]["rho"],
                REF["augmented_scvi_ord"]["rho"],
                ext["rho_overall"]]
        x = np.arange(len(methods))
        w = 0.35
        c_qwk = [MASLD_PAL["vanilla_ord"], MASLD_PAL["augmented_ord"],
                 MASLD_PAL["scanvi_aug"]]
        ax.bar(x - w/2, qwks, w, color=c_qwk, alpha=0.9, label="Andrews CV QWK")
        ax.bar(x + w/2, rhos, w,
               color=[MASLD_PAL["vanilla_ord"], MASLD_PAL["augmented_ord"],
                      MASLD_PAL["scanvi_aug"]],
               hatch="//", alpha=0.6,
               label=r"External $\rho$ vs disease_stage")
        for i, (q, r) in enumerate(zip(qwks, rhos)):
            if not np.isnan(q):
                ax.text(i - w/2, q + 0.02, f"{q:.3f}", ha="center", fontsize=8)
            ax.text(i + w/2, r + 0.02, f"{r:.3f}", ha="center", fontsize=8)
        ax.axhline(0.639, ls="--", lw=0.7, color="black", alpha=0.5,
                    label=r"target $\rho \geq 0.639$")
        ax.set_xticks(x)
        ax.set_xticklabels(methods, fontsize=8)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("Score")
        ax.set_title("(a) Head-to-head: QWK and external Spearman")
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(axis="y", alpha=0.25)

        # Panel B: SH breakdown (predicted F-stage distribution)
        ax = axes[0, 1]
        stages = ["F0", "F1", "F2", "F3", "F4"]
        x = np.arange(5)
        w = 0.27
        ax.bar(x - w, REF["vanilla_scvi_ord"]["sh_frac"], w,
               color=MASLD_PAL["vanilla_ord"],
               label="vanilla scVI cell-frac")
        ax.bar(x, REF["augmented_scvi_ord"]["sh_frac"], w,
               color=MASLD_PAL["augmented_ord"],
               label="augmented scVI ordinal (343m)")
        ax.bar(x + w, ext["sh_frac"], w,
               color=MASLD_PAL["scanvi_aug"],
               label="scANVI augmented (this)")
        ax.set_xticks(x)
        ax.set_xticklabels(stages)
        ax.set_ylim(0, 0.6)
        ax.set_ylabel("Fraction of SH donors")
        ax.set_title(f"(b) Steatohepatitis donor F-stage (n={ext['n_sh']})")
        ax.legend(fontsize=8)
        ax.grid(axis="y", alpha=0.25)

        # Panel C: per-dataset external rho
        ax = axes[1, 0]
        perds = ext["per_dataset"].dropna(subset=["rho"]).sort_values("rho")
        if len(perds) > 0:
            colors = [MASLD_PAL["scanvi_aug"] if r >= 0.5
                      else "#C44536" for r in perds["rho"].values]
            ax.barh(perds["dataset"].values, perds["rho"].values,
                    color=colors)
            for i, (rho, n) in enumerate(zip(perds["rho"].values,
                                              perds["n"].values)):
                ax.text(rho + 0.01, i, f"{rho:.2f} (n={n})",
                         va="center", fontsize=8)
            ax.axvline(0.5, ls="--", color="black", lw=0.7, alpha=0.6,
                       label="target rho >= 0.5")
            ax.axvline(0, ls="-", color="gray", lw=0.5)
            ax.set_xlim(-0.3, 1.15)
            ax.set_xlabel(r"Spearman $\rho$")
            ax.set_title("(c) Per-dataset external rho")
            ax.legend(fontsize=8, loc="lower right")
            ax.grid(axis="x", alpha=0.25)
        else:
            ax.text(0.5, 0.5, "no per-dataset data",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()

        # Panel D: Healthy->F4 hallucinations
        ax = axes[1, 1]
        labels_d = ["vanilla scVI\nordinal (343b)",
                    "augmented scVI\nordinal (343m)",
                    "scANVI\naugmented (this)"]
        vals = [REF["vanilla_scvi_ord"]["h2f4"],
                REF["augmented_scvi_ord"]["h2f4"],
                ext["n_healthy_f4"]]
        bars = ax.bar(np.arange(len(labels_d)), vals,
                      color=[MASLD_PAL["vanilla_ord"],
                             MASLD_PAL["augmented_ord"],
                             MASLD_PAL["scanvi_aug"]])
        for i, v in enumerate(vals):
            ax.text(i, v + 0.5, str(v), ha="center", fontsize=10,
                    fontweight="bold")
        ax.set_xticks(np.arange(len(labels_d)))
        ax.set_xticklabels(labels_d, fontsize=8)
        ax.set_ylabel("Healthy donors predicted F4")
        ax.set_title("(d) Hallucination check: clean-healthy -> F4 errors")
        ax.axhline(2, ls=":", color="black", lw=0.7, alpha=0.5,
                    label="target <= 2")
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(axis="y", alpha=0.25)
        ax.set_ylim(0, max(vals) * 1.3 + 1)

        fig.suptitle(
            f"scANVI with augmented F-stage labels  "
            f"(n_labelled_cells={n_train_cells:,})",
            fontsize=12, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
    log(f"  wrote {OUT_PDF}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subsample", type=int, default=0,
                    help="Subsample atlas to N cells (stratified by dataset).")
    ap.add_argument("--retrain-scvi", action="store_true",
                    help="Retrain scVI from scratch (default: load persisted).")
    ap.add_argument("--skip-cv", action="store_true",
                    help="Skip CV (faster: just train once and predict).")
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
    fstage["F_stage_documented"] = pd.to_numeric(
        fstage["F_stage_documented"], errors="coerce")
    meta = pd.read_csv(META_EXT, sep="\t")
    dubious = pd.read_csv(DUBIOUS, sep="\t")
    dubious_samples = set(dubious["sample"].astype(str).values)
    log(f"loaded: {fstage['F_stage_documented'].notna().sum()} documented "
        f"Andrews; {len(meta)} total donors in meta; "
        f"{len(dubious_samples)} dubious-healthy donors")

    # Cell-level augmented labels
    cell_labels = broadcast_augmented_labels(adata, fstage, meta, dubious_samples)
    n_labelled = int((cell_labels != UNLABELED).sum())
    adata.obs["fstage_aug_label"] = pd.Categorical(
        cell_labels, categories=F_STAGE_STR + [UNLABELED])

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
                n_latent=SCVI_N_LATENT, n_layers=SCVI_N_LAYERS,
                n_hidden=SCVI_N_HIDDEN,
                gene_likelihood=SCVI_GENE_LIKELIHOOD, dropout_rate=SCVI_DROPOUT,
            )
            scvi_model.train(
                max_epochs=150, early_stopping=True, early_stopping_patience=10,
                train_size=0.9, batch_size=512, enable_progress_bar=False,
            )
            scvi_model.save(str(MODEL_OUT / "scvi_full"), overwrite=True)
    else:
        log("training scVI from scratch")
        scvi_model = scvi.model.SCVI(
            adata,
            n_latent=SCVI_N_LATENT, n_layers=SCVI_N_LAYERS,
            n_hidden=SCVI_N_HIDDEN,
            gene_likelihood=SCVI_GENE_LIKELIHOOD, dropout_rate=SCVI_DROPOUT,
        )
        scvi_model.train(
            max_epochs=150, early_stopping=True, early_stopping_patience=10,
            train_size=0.9, batch_size=512, enable_progress_bar=False,
        )
        scvi_model.save(str(MODEL_OUT / "scvi_full"), overwrite=True)

    # ── K-fold CV (within Andrews; anchors stay in training) ──────────────
    cv_result = None
    if not args.skip_cv:
        log(f"===== {args.cv_folds}-fold CV (Andrews held out; "
            f"clean-healthy and cirrhosis anchors stay in training) =====")
        try:
            cv_result = fold_cv(
                adata, fstage, meta, dubious_samples,
                cell_labels, scvi_model,
                n_splits=args.cv_folds, seed=42,
            )
        except Exception as e:
            log(f"CV FAILED: {e}")
            import traceback; traceback.print_exc()
            cv_result = None

    # ── Final scANVI: train on FULL augmented labels and predict every cell ──
    log("===== final scANVI on full augmented labels =====")
    adata.obs["fstage_aug_label"] = pd.Categorical(
        cell_labels, categories=F_STAGE_STR + [UNLABELED])
    scanvi_full = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        labels_key="fstage_aug_label",
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
    log(f"  saved scANVI augmented -> {MODEL_OUT / 'scanvi_full'}")

    # ── Predict for every hepatocyte ──────────────────────────────────────
    log("predicting F-stage for all hepatocytes")
    argmax_int, proba, classes = predict_with_proba(scanvi_full, adata)
    per_cell = pd.DataFrame({
        "cell_id": adata.obs_names.values.astype(str),
        "sample":  adata.obs["sample"].astype(str).values,
        "dataset": adata.obs["dataset"].astype(str).values,
        "F_stage_scanvi_aug": argmax_int,
    })
    for j, c in enumerate(classes):
        per_cell[f"P_F{c}"] = proba[:, j].astype(np.float32)
    per_cell["F_stage"] = per_cell["F_stage_scanvi_aug"]  # alias per spec
    # write per-cell (compressed)
    out_percell_cols = ["cell_id", "sample", "dataset", "F_stage"] + \
                       [f"P_F{c}" for c in classes]
    per_cell[out_percell_cols].to_csv(
        OUT_PERCELL, sep="\t", index=False,
        compression="gzip", float_format="%.5f")
    log(f"  wrote {OUT_PERCELL} ({len(per_cell):,} rows)")

    # ── Donor-level aggregation ────────────────────────────────────────────
    p_cols = [f"P_F{c}" for c in classes]
    donor_agg = aggregate_donor(per_cell, p_cols)
    donor_agg.to_csv(OUT_DONOR_TSV, sep="\t", index=False, float_format="%.6f")
    log(f"  wrote {OUT_DONOR_TSV} ({len(donor_agg):,} donors)")

    # ── External validation ──────────────────────────────────────────────
    ext = external_validation(donor_agg, meta)

    # ── Method metrics + comparison append ───────────────────────────────
    qwk = cv_result["qwk"] if cv_result else float("nan")
    acc = cv_result["accuracy"] if cv_result else float("nan")
    n_train = n_labelled
    off_by_one = ""
    if cv_result is not None:
        cv_preds = cv_result["cv"]
        diff = (cv_preds["F_stage_documented"].astype(int) -
                cv_preds["F_stage_scanvi_augmented_argmax"].astype(int)).abs()
        off_by_one = round(float((diff <= 1).mean()), 4)
    method_row = dict(
        method="scanvi_augmented",
        n_train_cells=n_train,
        n_andrews_donors=int(fstage["F_stage_documented"].notna().sum()),
        qwk_andrews_loocv=round(qwk, 4) if cv_result else "",
        accuracy_andrews_loocv=round(acc, 4) if cv_result else "",
        off_by_one_loocv=off_by_one,
        external_rho=round(ext["rho_overall"], 4),
        external_rho_p="{:.3e}".format(ext["p_overall"]),
        n_healthy_F4_misassign=ext["n_healthy_f4"],
        cv_strategy=f"{args.cv_folds}-fold stratified Andrews" if cv_result else "skip",
        n_cells=int(adata.n_obs),
        sh_pred_F0_frac=ext["sh_frac"][0],
        sh_pred_F1_frac=ext["sh_frac"][1],
        sh_pred_F2_frac=ext["sh_frac"][2],
        sh_pred_F3_frac=ext["sh_frac"][3],
        sh_pred_F4_frac=ext["sh_frac"][4],
        n_sh_pred=ext["n_sh"],
    )
    pd.DataFrame([method_row]).to_csv(OUT_METRICS, sep="\t", index=False)
    log(f"  wrote {OUT_METRICS}")

    # Append to fstage_method_comparison.tsv (preserve schema)
    if OUT_COMPARE.exists():
        comp = pd.read_csv(OUT_COMPARE, sep="\t")
    else:
        comp = pd.DataFrame(columns=["method", "n_train", "qwk_loocv",
                                     "accuracy_loocv", "off_by_one_loocv"])
    comp = comp[comp["method"] != "scanvi_augmented"].copy()
    new_row = {
        "method": "scanvi_augmented",
        "n_train": int(fstage["F_stage_documented"].notna().sum()),
        "qwk_loocv": round(qwk, 4) if cv_result else "",
        "accuracy_loocv": round(acc, 4) if cv_result else "",
        "off_by_one_loocv": off_by_one,
    }
    comp = pd.concat([comp, pd.DataFrame([new_row])], ignore_index=True)
    comp.to_csv(OUT_COMPARE, sep="\t", index=False)
    log(f"  appended scanvi_augmented row -> {OUT_COMPARE}")

    # ── External validation summary ──────────────────────────────────────
    with open(OUT_SUMMARY, "w") as f:
        f.write("# scANVI AUGMENTED F-stage external validation\n")
        f.write(f"# cells: {adata.n_obs:,}; labelled cells: {n_labelled:,}; "
                f"Andrews documented donors: "
                f"{int(fstage['F_stage_documented'].notna().sum())}\n")
        f.write(f"# CV: {args.cv_folds}-fold stratified across Andrews "
                f"documented donors (anchors stay in training)\n")
        if cv_result is not None:
            f.write(f"\nCV QWK    = {cv_result['qwk']:.3f}\n")
            f.write(f"CV accuracy = {cv_result['accuracy']:.3f}\n")
            f.write(f"CV off-by-one tolerance = {off_by_one}\n")
            f.write("\nConfusion matrix (rows=true, cols=pred, F0..F4):\n")
            f.write(str(np.array(cv_result["cm"])) + "\n")
        f.write(f"\nExternal Spearman rho(F_stage_scanvi_aug, "
                f"disease_stage_numeric) = {ext['rho_overall']:.3f} "
                f"(p={ext['p_overall']:.2e})\n")
        f.write(f"Healthy donors predicted as F4: {ext['n_healthy_f4']}\n")
        f.write("\nPer-dataset Spearman rho:\n")
        f.write(ext["per_dataset"].to_string(index=False) + "\n")
        f.write("\nPredicted F-stage x disease_stage_coarse:\n")
        f.write(ext["crosstab"].to_string() + "\n")
        f.write(f"\nSteatohepatitis predicted F-stage (n={ext['n_sh']}):\n")
        for k, v in enumerate(ext["sh_frac"]):
            f.write(f"  F{k}: {v*100:5.1f}%  (n={ext['sh_counts'][k]})\n")
        f.write("\nHead-to-head:\n")
        f.write(f"  vanilla scVI ordinal (343b):   QWK=0.7422  rho=0.391  H->F4=27\n")
        f.write(f"  augmented scVI ordinal (343m): QWK=0.7200  rho=0.639  H->F4=2\n")
        f.write(f"  scANVI augmented (this):       "
                f"QWK={qwk if cv_result else float('nan'):.4f}  "
                f"rho={ext['rho_overall']:.4f}  "
                f"H->F4={ext['n_healthy_f4']}\n")
        rho_target_met = ext["rho_overall"] >= 0.639
        f4_target_met = ext["n_healthy_f4"] <= 2
        per_ds = ext["per_dataset"]
        target_datasets = {"GSE136103", "GSE174748", "GSE244832"}
        per_target = per_ds[per_ds["dataset"].isin(target_datasets)]
        per_target_ok = (per_target.dropna(subset=["rho"])["rho"] >= 0.5).all() \
            if len(per_target) else False
        f.write("\nSuccess criteria:\n")
        f.write(f"  external rho >= 0.639:                 "
                f"{'PASS' if rho_target_met else 'FAIL'}\n")
        f.write(f"  Healthy->F4 errors <= 2:               "
                f"{'PASS' if f4_target_met else 'FAIL'}\n")
        f.write(f"  per-dataset rho >= 0.5 (target trio):  "
                f"{'PASS' if per_target_ok else 'FAIL'}\n")
        f.write(f"\nTarget datasets per-dataset rho:\n")
        if len(per_target):
            f.write(per_target.to_string(index=False) + "\n")
    log(f"  wrote {OUT_SUMMARY}")

    # ── Figure ───────────────────────────────────────────────────────────
    make_figure(cv_result, ext, donor_agg, meta, n_labelled)

    log("DONE")


if __name__ == "__main__":
    main()
