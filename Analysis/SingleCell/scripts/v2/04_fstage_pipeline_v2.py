#!/usr/bin/env python
"""
04_fstage_pipeline_v2.py

Phase 0.5 v2 equivalent of the 343b / 343m / 343e / 343o F-stage suite.

Trains an ordinal-logistic F-stage classifier on donor-level scVI v2 pseudobulks
and emits four artifacts under
  Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2/

  1. donor_fstage_scvi_v2_predicted.tsv      (343b equivalent)
  2. donor_fstage_augmented_v2_predicted.tsv (343m equivalent; with anchors)
  3. percell_fstage_predictions_v2.tsv.gz +  (343e equivalent)
     donor_fstage_cellfrac_v2.tsv
  4. heldout_cohort_metrics_v2.tsv +         (343o equivalent)
     heldout_cohort_fstage_replication_v2.tsv

A consolidated fstage_method_comparison_v2.tsv table records v2 metrics with
the `v2` tag for head-to-head against the v1 numbers.

This script reads ONLY:
  - Agent 3's hepatocyte_atlas_v2_annotated.h5ad (v2 hepatocyte atlas)
  - V1 donor_fstage_documented.tsv (read-only documented labels)
  - V1 dubious_healthy_donors.tsv (read-only — used to gate clean-healthy anchors)
  - V1 donor_metadata.tsv (read-only — provides disease_stage_coarse roster)

All outputs are written under results_gpu_v2_phase05; the script NEVER touches
v1 outputs.
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

from sklearn.metrics import (
    cohen_kappa_score,
    confusion_matrix,
    accuracy_score,
    roc_auc_score,
)
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from scipy import stats

try:
    import mord

    HAVE_MORD = True
except ImportError:
    HAVE_MORD = False

# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------
PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)

# v2 inputs (Agent 3 product) — READ-WRITE on v2 outputs
HEP_H5AD_V2 = (
    PROJECT_ROOT
    / "Analysis/SingleCell/results_gpu_v2_phase05/atlas/hepatocyte_atlas_v2_annotated.h5ad"
)

# v1 reference (read-only)
V1_STAGE_DIR = (
    PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
)
FSTAGE_DOC_V1 = V1_STAGE_DIR / "donor_fstage_documented.tsv"
DUBIOUS_V1 = V1_STAGE_DIR / "dubious_healthy_donors.tsv"
DONOR_META_V1 = (
    PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv"
)

# v2 outputs
OUT_DIR = (
    PROJECT_ROOT
    / "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2"
)
OUT_SCVI = OUT_DIR / "donor_fstage_scvi_v2_predicted.tsv"
OUT_AUGMENTED = OUT_DIR / "donor_fstage_augmented_v2_predicted.tsv"
OUT_PERCELL = OUT_DIR / "percell_fstage_predictions_v2.tsv.gz"
OUT_DONOR_FRAC = OUT_DIR / "donor_fstage_cellfrac_v2.tsv"
OUT_HEP_PB = OUT_DIR / "donor_hep_scvi_mean_v2.tsv"
OUT_HELDOUT = OUT_DIR / "heldout_cohort_metrics_v2.tsv"
OUT_HELDOUT_PRED = OUT_DIR / "heldout_cohort_fstage_replication_v2.tsv"
OUT_METHOD_CMP = OUT_DIR / "fstage_method_comparison_v2.tsv"

QWK_HALT_THRESHOLD = 0.10  # below this, latent space lacks F-stage signal
AUROC_GENERALIZES_THRESHOLD = 0.65

# Pre-registered cohorts to evaluate held-out replication. We discover the
# actual evaluable cohort list at runtime since v2 may have fewer cohorts
# with sufficient donors than v1.
HELDOUT_MIN_DONORS = 3


# ======================================================================
# Helpers
# ======================================================================
def detect_latent_dim_and_obsm(a) -> tuple[str, int, list[str]]:
    """Return (obsm_key, latent_dim, latent_cols) for the v2 scVI latent.

    Tries common keys: X_scVI, X_scanvi, X_scvi_v2. Halts if none found.
    """
    candidate_keys = ["X_scVI", "X_scanvi", "X_scvi_v2", "X_scanvi_v2"]
    for k in candidate_keys:
        if k in a.obsm:
            Z = np.asarray(a.obsm[k])
            dim = Z.shape[1]
            cols = [f"z{i}" for i in range(dim)]
            print(f"[detect] using obsm['{k}'] (n_cells={Z.shape[0]}, dim={dim})")
            return k, dim, cols
    sys.exit(
        f"FATAL: no recognized scVI latent in obsm; keys={list(a.obsm.keys())}"
    )


def make_ordinal():
    if HAVE_MORD:
        return Pipeline(
            [
                ("scale", StandardScaler()),
                ("clf", mord.LogisticIT(alpha=1.0)),
            ]
        )
    return Pipeline(
        [
            ("scale", StandardScaler()),
            (
                "clf",
                LogisticRegression(
                    solver="lbfgs",
                    C=1.0,
                    max_iter=5000,
                    class_weight="balanced",
                ),
            ),
        ]
    )


def predict_argmax(clf, X):
    if HAVE_MORD and isinstance(clf.named_steps["clf"], mord.LogisticIT):
        return clf.predict(X).astype(int)
    proba = clf.predict_proba(X)
    return clf.named_steps["clf"].classes_[proba.argmax(axis=1)].astype(int)


def predict_proba_full(clf, X) -> np.ndarray:
    proba = clf.predict_proba(X)
    classes = clf.named_steps["clf"].classes_
    full = np.zeros((proba.shape[0], 5), dtype=float)
    for j, c in enumerate(classes):
        full[:, int(c)] = proba[:, j]
    rs = full.sum(axis=1, keepdims=True)
    rs[rs == 0] = 1.0
    return full / rs


def loocv_qwk(X: np.ndarray, y: np.ndarray) -> tuple[float, np.ndarray]:
    n = len(y)
    yhat = np.empty(n, dtype=int)
    for i in range(n):
        m = np.ones(n, dtype=bool)
        m[i] = False
        clf = make_ordinal()
        clf.fit(X[m], y[m])
        yhat[i] = predict_argmax(clf, X[i : i + 1])[0]
    qwk = cohen_kappa_score(
        y, yhat, weights="quadratic", labels=list(range(5))
    )
    return float(qwk), yhat


# ======================================================================
# Data assembly
# ======================================================================
def build_donor_pb(
    a: ad.AnnData, obsm_key: str, latent_cols: list[str]
) -> pd.DataFrame:
    """Donor-level mean of scVI latent over hepatocytes."""
    print(f"[pb] computing donor pseudobulk from obsm['{obsm_key}']")
    Z = np.asarray(a.obsm[obsm_key])
    samples = a.obs["sample"].astype(str).to_numpy()
    datasets = a.obs["dataset"].astype(str).to_numpy()

    df = pd.DataFrame(Z, columns=latent_cols)
    df["sample"] = samples
    df["dataset"] = datasets
    pb = df.groupby(["sample", "dataset"], observed=True)[
        latent_cols
    ].mean().reset_index()
    n_hep = (
        df.groupby(["sample", "dataset"], observed=True)
        .size()
        .rename("n_hepatocytes_for_latent")
        .reset_index()
    )
    pb = pb.merge(n_hep, on=["sample", "dataset"])
    print(f"[pb] {len(pb)} (sample, dataset) rows")
    return pb


def load_meta(pb: pd.DataFrame, latent_cols: list[str]) -> pd.DataFrame:
    """Join v1 donor_metadata + documented F-stage + dubious flag into pb."""
    meta = pd.read_csv(DONOR_META_V1, sep="\t")
    doc = pd.read_csv(FSTAGE_DOC_V1, sep="\t")
    doc["F_stage_documented"] = pd.to_numeric(
        doc["F_stage_documented"], errors="coerce"
    )
    dub_path = DUBIOUS_V1
    if dub_path.exists():
        dub = pd.read_csv(dub_path, sep="\t")
        dubious_samples = set(dub["sample"].astype(str).values)
    else:
        print(f"[warn] no dubious file at {dub_path}; treating all Healthy as clean")
        dubious_samples = set()

    keep_meta = ["sample", "dataset", "disease_stage_coarse"]
    if "disease_stage_numeric" in meta.columns:
        keep_meta.append("disease_stage_numeric")
    df = meta[keep_meta].copy()
    df = df.merge(pb, on=["sample", "dataset"], how="left")
    df = df.merge(
        doc[["sample", "dataset", "F_stage_documented"]],
        on=["sample", "dataset"],
        how="left",
    )
    df["dubious"] = df["sample"].astype(str).isin(dubious_samples)
    if "disease_stage_numeric" not in df.columns:
        stage_map = {
            "Healthy": 0,
            "Steatosis": 1,
            "Steatohepatitis": 2,
            "Cirrhosis": 3,
        }
        df["disease_stage_numeric"] = df["disease_stage_coarse"].map(stage_map)
    return df


def assemble_anchors(df: pd.DataFrame, latent_cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    out["y_train"] = np.nan
    out["origin"] = "none"

    mask_doc = out["F_stage_documented"].notna()
    out.loc[mask_doc, "y_train"] = out.loc[mask_doc, "F_stage_documented"]
    out.loc[mask_doc, "origin"] = "documented"

    mask_ch = (
        (out["disease_stage_coarse"] == "Healthy")
        & (~out["dubious"])
        & (~mask_doc)
    )
    out.loc[mask_ch, "y_train"] = 0
    out.loc[mask_ch, "origin"] = "clean_healthy_anchor"

    mask_cirr = (out["disease_stage_coarse"] == "Cirrhosis") & (~mask_doc)
    out.loc[mask_cirr, "y_train"] = 4
    out.loc[mask_cirr, "origin"] = "cirrhosis_anchor"

    return out


# ======================================================================
# (1) 343b equivalent: documented-only training
# ======================================================================
def run_scvi_only(
    df: pd.DataFrame, latent_cols: list[str]
) -> dict:
    print("\n========== (343b equiv): documented-only training ==========")
    train = df[df["F_stage_documented"].notna() & df[latent_cols].notna().all(axis=1)].copy()
    train["F_stage_documented"] = train["F_stage_documented"].astype(int)
    if len(train) < 20:
        sys.exit(f"FATAL: only {len(train)} documented donors with latents")

    X = train[latent_cols].to_numpy()
    y = train["F_stage_documented"].to_numpy().astype(int)
    print(f"[train] n={len(y)} class counts: {dict(zip(*np.unique(y, return_counts=True)))}")
    qwk, yhat = loocv_qwk(X, y)
    print(f"[loocv] QWK = {qwk:.4f}")
    print("[loocv] CM (rows=true F0..F4, cols=pred):")
    print(confusion_matrix(y, yhat, labels=list(range(5))))

    if qwk < QWK_HALT_THRESHOLD:
        print(
            f"[warn] QWK {qwk:.4f} < {QWK_HALT_THRESHOLD}; writing anyway "
            "so downstream can see degradation."
        )

    # Refit on full training and predict everyone with latents
    clf = make_ordinal()
    clf.fit(X, y)
    pred_mask = df[latent_cols].notna().all(axis=1)
    X_all = df.loc[pred_mask, latent_cols].to_numpy()
    pred = predict_argmax(clf, X_all)
    proba = predict_proba_full(clf, X_all)

    out = df[["sample", "dataset"]].copy()
    out["F_stage_predicted_argmax"] = np.nan
    out.loc[pred_mask, "F_stage_predicted_argmax"] = pred.astype(float)
    for k in range(5):
        col = f"P_F{k}"
        out[col] = np.nan
        out.loc[pred_mask, col] = proba[:, k].round(6)
    out["classifier_qwk_loocv"] = round(qwk, 4)
    out["n_hepatocytes"] = df["n_hepatocytes_for_latent"]

    out.to_csv(OUT_SCVI, sep="\t", index=False)
    print(f"[output] {OUT_SCVI}")

    return {
        "method": "scvi_v2_ordinal",
        "n_train": int(len(y)),
        "qwk_loocv": float(qwk),
        "version": "v2",
    }


# ======================================================================
# (2) 343m equivalent: augmented training (anchors)
# ======================================================================
def run_augmented(
    df_anch: pd.DataFrame, latent_cols: list[str]
) -> tuple[dict, pd.Series]:
    print("\n========== (343m equiv): augmented training ==========")
    train = df_anch[
        df_anch["y_train"].notna()
        & df_anch[latent_cols].notna().all(axis=1)
    ].copy()
    train["y_train"] = train["y_train"].astype(int)
    print(f"[train] n_total={len(train)} composition:")
    print(train.groupby(["origin", "y_train"]).size().unstack(fill_value=0))

    X = train[latent_cols].to_numpy()
    y = train["y_train"].to_numpy().astype(int)
    origin = train["origin"].to_numpy()

    qwk_full, yhat = loocv_qwk(X, y)
    print(f"[full-LOOCV] QWK = {qwk_full:.4f}")

    # Andrews-only LOOCV (apples-to-apples vs 343b)
    andrews_idx = np.where(origin == "documented")[0]
    if len(andrews_idx) > 0:
        yhat_a = np.full(len(y), -1, dtype=int)
        for i in andrews_idx:
            m = np.ones(len(y), dtype=bool)
            m[i] = False
            clf = make_ordinal()
            clf.fit(X[m], y[m])
            yhat_a[i] = predict_argmax(clf, X[i : i + 1])[0]
        qwk_andrews = cohen_kappa_score(
            y[andrews_idx],
            yhat_a[andrews_idx],
            weights="quadratic",
            labels=list(range(5)),
        )
        print(
            f"[andrews-LOOCV] QWK = {qwk_andrews:.4f} (over {len(andrews_idx)} docs)"
        )
    else:
        qwk_andrews = float("nan")

    # Healthy->F4 hallucinations
    ch_idx = np.where(origin == "clean_healthy_anchor")[0]
    n_h_to_F4 = 0
    if len(ch_idx) > 0:
        yhat_ch = np.full(len(y), -1, dtype=int)
        for i in ch_idx:
            m = np.ones(len(y), dtype=bool)
            m[i] = False
            clf = make_ordinal()
            clf.fit(X[m], y[m])
            yhat_ch[i] = predict_argmax(clf, X[i : i + 1])[0]
        n_h_to_F4 = int((yhat_ch[ch_idx] == 4).sum())
        print(
            f"[hallu] clean-healthy LOOCV: {n_h_to_F4}/{len(ch_idx)} predicted F4"
        )

    # Refit on full anchors set, predict everyone with latents
    clf = make_ordinal()
    clf.fit(X, y)
    pred_mask = df_anch[latent_cols].notna().all(axis=1)
    X_all = df_anch.loc[pred_mask, latent_cols].to_numpy()
    pred = predict_argmax(clf, X_all)
    proba = predict_proba_full(clf, X_all)

    out = df_anch[
        ["sample", "dataset", "disease_stage_coarse", "disease_stage_numeric"]
    ].copy()
    out["F_stage_documented"] = df_anch["F_stage_documented"]
    out["origin"] = df_anch["origin"]
    out["F_stage_pred_augmented_scvi_v2_ordinal"] = np.nan
    out.loc[
        pred_mask, "F_stage_pred_augmented_scvi_v2_ordinal"
    ] = pred.astype(float)
    for k in range(5):
        col = f"P_F{k}_augmented_scvi_v2_ordinal"
        out[col] = np.nan
        out.loc[pred_mask, col] = proba[:, k].round(6)

    # External Spearman vs disease_stage_numeric
    valid = out.dropna(
        subset=[
            "F_stage_pred_augmented_scvi_v2_ordinal",
            "disease_stage_numeric",
        ]
    )
    rho, rho_p = (
        stats.spearmanr(
            valid["F_stage_pred_augmented_scvi_v2_ordinal"],
            valid["disease_stage_numeric"],
        )
        if len(valid) >= 3
        else (np.nan, np.nan)
    )
    print(f"[external] Spearman rho vs disease_stage_numeric: {rho:.4f} (p={rho_p:.3e}, n={len(valid)})")

    out.to_csv(OUT_AUGMENTED, sep="\t", index=False)
    print(f"[output] {OUT_AUGMENTED}")

    augmented_pred_series = out.set_index(["sample", "dataset"])[
        "F_stage_pred_augmented_scvi_v2_ordinal"
    ]

    return (
        {
            "method": "augmented_scvi_v2_ordinal",
            "n_train": int(len(y)),
            "qwk_full_loocv": float(qwk_full),
            "qwk_andrews_loocv": float(qwk_andrews),
            "external_rho": float(rho) if not np.isnan(rho) else None,
            "external_rho_p": float(rho_p) if not np.isnan(rho_p) else None,
            "n_external": int(len(valid)),
            "n_healthy_to_F4": int(n_h_to_F4),
            "version": "v2",
        },
        augmented_pred_series,
    )


# ======================================================================
# (3) 343e equivalent: per-cell breakdown
# ======================================================================
def run_percell(
    a: ad.AnnData,
    obsm_key: str,
    latent_cols: list[str],
    df: pd.DataFrame,
) -> dict:
    print("\n========== (343e equiv): per-cell F-stage breakdown ==========")
    # Train LR on Andrews-only pseudobulks (same recipe as 343e)
    train = df[
        df["F_stage_documented"].notna() & df[latent_cols].notna().all(axis=1)
    ].copy()
    train["F_stage_documented"] = train["F_stage_documented"].astype(int)
    X = train[latent_cols].to_numpy(dtype=np.float64)
    y = train["F_stage_documented"].to_numpy()
    scaler = StandardScaler().fit(X)
    clf = LogisticRegression(
        solver="lbfgs", C=1.0, max_iter=5000, class_weight="balanced"
    )
    clf.fit(scaler.transform(X), y)
    print(f"[train] n={len(y)}, classes={clf.classes_.tolist()}")

    # Per-cell prediction in chunks
    Z = np.asarray(a.obsm[obsm_key], dtype=np.float32)
    samples = a.obs["sample"].astype(str).to_numpy()
    datasets = a.obs["dataset"].astype(str).to_numpy()
    n = Z.shape[0]
    chunk = 200_000
    argmax = np.empty(n, dtype=np.int8)
    proba_full = np.empty((n, 5), dtype=np.float32)
    for i in range(0, n, chunk):
        j = min(i + chunk, n)
        Xc = scaler.transform(Z[i:j].astype(np.float64))
        p = clf.predict_proba(Xc)
        full = np.zeros((p.shape[0], 5), dtype=np.float32)
        for jj, c in enumerate(clf.classes_):
            full[:, int(c)] = p[:, jj]
        rs = full.sum(axis=1, keepdims=True)
        rs[rs == 0] = 1.0
        full = full / rs
        argmax[i:j] = full.argmax(axis=1).astype(np.int8)
        proba_full[i:j] = full
        print(f"[predict] cells {i}-{j}/{n}", flush=True)

    percell = pd.DataFrame(
        {
            "sample": samples,
            "dataset": datasets,
            "F_stage_argmax": argmax.astype(int),
            "P_F0": proba_full[:, 0].round(5),
            "P_F1": proba_full[:, 1].round(5),
            "P_F2": proba_full[:, 2].round(5),
            "P_F3": proba_full[:, 3].round(5),
            "P_F4": proba_full[:, 4].round(5),
        }
    )
    percell.to_csv(OUT_PERCELL, sep="\t", index=False, compression="gzip")
    print(f"[output] {OUT_PERCELL} ({len(percell)} rows)")

    # Donor cellfrac
    counts = (
        percell.groupby(["sample", "dataset", "F_stage_argmax"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    for k in range(5):
        if k not in counts.columns:
            counts[k] = 0
    counts = counts[[0, 1, 2, 3, 4]]
    counts.columns = [f"F{k}_count" for k in range(5)]
    counts["n_cells"] = counts.sum(axis=1)
    fracs = counts[[f"F{k}_count" for k in range(5)]].div(
        counts["n_cells"], axis=0
    )
    fracs.columns = [f"F{k}_frac" for k in range(5)]
    donor_frac = pd.concat([counts, fracs], axis=1).reset_index()
    donor_frac["dominant_fstage"] = fracs.values.argmax(axis=1)

    def _entropy(v):
        v = np.clip(v, 1e-12, 1.0)
        v = v / v.sum()
        return float(-np.sum(v * np.log2(v)))

    donor_frac["entropy"] = [_entropy(v) for v in fracs.values]
    donor_frac.to_csv(OUT_DONOR_FRAC, sep="\t", index=False)
    print(f"[output] {OUT_DONOR_FRAC}")

    return {
        "method": "percell_v2",
        "n_train": int(len(y)),
        "n_cells_scored": int(n),
        "version": "v2",
    }


# ======================================================================
# (4) 343o equivalent: held-out cohort replication
# ======================================================================
def run_heldout(df_anch: pd.DataFrame, latent_cols: list[str]) -> list[dict]:
    print("\n========== (343o equiv): held-out cohort replication ==========")
    cohort_counts = (
        df_anch[df_anch["y_train"].notna()].groupby("dataset").size()
    )
    eval_datasets = [
        ds for ds, n in cohort_counts.items() if n >= HELDOUT_MIN_DONORS
    ]
    print(f"[loo] {len(eval_datasets)} evaluable datasets: {eval_datasets}")

    out_rows = []
    pred_rows = []
    for ds in eval_datasets:
        train = df_anch[
            df_anch["y_train"].notna()
            & df_anch[latent_cols].notna().all(axis=1)
            & (df_anch["dataset"] != ds)
        ].copy()
        train["y_train"] = train["y_train"].astype(int)
        hold = df_anch[
            (df_anch["dataset"] == ds)
            & df_anch[latent_cols].notna().all(axis=1)
        ].copy()
        if len(train) < 10 or len(hold) == 0:
            print(f"[skip:{ds}] insufficient data")
            continue
        X_tr = train[latent_cols].to_numpy()
        y_tr = train["y_train"].to_numpy().astype(int)
        clf = make_ordinal()
        clf.fit(X_tr, y_tr)
        X_ho = hold[latent_cols].to_numpy()
        proba = predict_proba_full(clf, X_ho)
        pred = proba.argmax(axis=1)

        # Per-donor pred frame
        hold = hold.assign(
            F_stage_pred=pred,
            P_F0=proba[:, 0],
            P_F1=proba[:, 1],
            P_F2=proba[:, 2],
            P_F3=proba[:, 3],
            P_F4=proba[:, 4],
        )
        hold["held_out_cohort"] = ds
        pred_cols = [
            "held_out_cohort",
            "sample",
            "dataset",
            "disease_stage_coarse",
            "disease_stage_numeric",
            "F_stage_documented",
            "F_stage_pred",
            "P_F0",
            "P_F1",
            "P_F2",
            "P_F3",
            "P_F4",
        ]
        pred_rows.append(hold[[c for c in pred_cols if c in hold.columns]])

        # Metrics
        y_binary = (hold["disease_stage_coarse"] != "Healthy").astype(int)
        score_diseased = hold["P_F2"] + hold["P_F3"] + hold["P_F4"]
        auroc = (
            roc_auc_score(y_binary, score_diseased)
            if y_binary.nunique() == 2
            else np.nan
        )
        valid_rho = hold.dropna(subset=["disease_stage_numeric"])
        if (
            len(valid_rho) >= 3
            and valid_rho["disease_stage_numeric"].nunique() >= 2
        ):
            rho, rho_p = stats.spearmanr(
                valid_rho["F_stage_pred"], valid_rho["disease_stage_numeric"]
            )
        else:
            rho, rho_p = np.nan, np.nan
        n_h = int((hold["disease_stage_coarse"] == "Healthy").sum())
        n_h_to_F4 = int(
            (
                (hold["disease_stage_coarse"] == "Healthy")
                & (hold["F_stage_pred"] == 4)
            ).sum()
        )

        doc_mask = hold["F_stage_documented"].notna()
        qwk_doc = np.nan
        if doc_mask.sum() >= 3:
            qwk_doc = cohen_kappa_score(
                hold.loc[doc_mask, "F_stage_documented"].astype(int),
                hold.loc[doc_mask, "F_stage_pred"].astype(int),
                weights="quadratic",
                labels=list(range(5)),
            )

        if not np.isnan(auroc):
            passes = bool(auroc > AUROC_GENERALIZES_THRESHOLD)
            pmetric = "AUROC"
        elif not np.isnan(qwk_doc):
            passes = bool(qwk_doc > 0.5)
            pmetric = "QWK_doc"
        else:
            passes = False
            pmetric = "undefined"

        print(
            f"[{ds}] n_train={len(train)} n_hold={len(hold)} "
            f"AUROC={auroc if not np.isnan(auroc) else 'NA'} "
            f"rho={rho if not np.isnan(rho) else 'NA'} "
            f"QWK_doc={qwk_doc if not np.isnan(qwk_doc) else 'NA'} "
            f"pass={passes}"
        )
        out_rows.append(
            {
                "held_out_cohort": ds,
                "n_train": int(len(train)),
                "n_train_cohorts": int(train["dataset"].nunique()),
                "n_holdout": int(len(hold)),
                "n_holdout_healthy": n_h,
                "auroc_healthy_vs_diseased": (
                    round(float(auroc), 4) if not np.isnan(auroc) else np.nan
                ),
                "spearman_rho": (
                    round(float(rho), 4) if not np.isnan(rho) else np.nan
                ),
                "qwk_vs_documented": (
                    round(float(qwk_doc), 4) if not np.isnan(qwk_doc) else np.nan
                ),
                "n_healthy_to_F4": n_h_to_F4,
                "pass_metric": pmetric,
                "passes_threshold": passes,
                "version": "v2",
            }
        )

    if out_rows:
        pd.DataFrame(out_rows).to_csv(OUT_HELDOUT, sep="\t", index=False)
        print(f"[output] {OUT_HELDOUT}")
    if pred_rows:
        pd.concat(pred_rows, ignore_index=True).to_csv(
            OUT_HELDOUT_PRED, sep="\t", index=False
        )
        print(f"[output] {OUT_HELDOUT_PRED}")

    return out_rows


# ======================================================================
# Main
# ======================================================================
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"[info] PROJECT_ROOT = {PROJECT_ROOT}")
    print(f"[info] HEP_H5AD_V2 = {HEP_H5AD_V2}")
    print(f"[info] OUT_DIR     = {OUT_DIR}")
    print(f"[info] HAVE_MORD={HAVE_MORD}")

    if not HEP_H5AD_V2.exists():
        sys.exit(f"FATAL: v2 hepatocyte atlas not found at {HEP_H5AD_V2}")

    t0 = time.time()
    a = ad.read_h5ad(HEP_H5AD_V2, backed="r")
    print(f"[load] {HEP_H5AD_V2.name}: n_cells={a.n_obs} n_genes={a.n_vars}")
    if "sample" not in a.obs.columns or "dataset" not in a.obs.columns:
        sys.exit(
            f"FATAL: obs missing 'sample' or 'dataset'. obs cols={list(a.obs.columns)}"
        )

    obsm_key, latent_dim, latent_cols = detect_latent_dim_and_obsm(a)

    pb = build_donor_pb(a, obsm_key, latent_cols)
    pb.to_csv(OUT_HEP_PB, sep="\t", index=False)
    print(f"[output] {OUT_HEP_PB} ({len(pb)} donors)")

    df = load_meta(pb, latent_cols)
    n_with_latent = int(df[latent_cols].notna().all(axis=1).sum())
    print(f"[meta] {n_with_latent}/{len(df)} donors have v2 latents")

    df_anch = assemble_anchors(df, latent_cols)
    n_anch = int(df_anch["y_train"].notna().sum())
    print(f"[anchors] {n_anch} anchored donors:")
    print(
        df_anch[df_anch["y_train"].notna()]
        .groupby(["origin", "y_train"])
        .size()
        .unstack(fill_value=0)
    )

    # Run the four arms
    method_rows = []
    m1 = run_scvi_only(df, latent_cols)
    method_rows.append(m1)
    m2, _ = run_augmented(df_anch, latent_cols)
    method_rows.append(m2)
    m3 = run_percell(a, obsm_key, latent_cols, df)
    method_rows.append(m3)
    heldout = run_heldout(df_anch, latent_cols)

    # Cleanup AnnData handle
    a.file.close()

    # Method comparison table
    cmp = pd.DataFrame(method_rows)
    cmp["latent_dim"] = latent_dim
    cmp.to_csv(OUT_METHOD_CMP, sep="\t", index=False)
    print(f"[output] {OUT_METHOD_CMP}")

    print(f"\n[done] 04_fstage_pipeline_v2 in {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
