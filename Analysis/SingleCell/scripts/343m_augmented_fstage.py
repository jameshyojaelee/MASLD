#!/usr/bin/env python
"""
343m_augmented_fstage.py

Augmented-training F-stage classifier.

Hypothesis
----------
Vanilla scVI F-stage classifier (343b) trains on 58 Andrews documented donors
only (LOOCV QWK = 0.742, external rho = 0.391). Failure mode: Andrews is the
only training cohort, so the classifier has no cross-cohort coverage and
hallucinates Healthy -> F4 (n=27 vanilla).

Cross-dataset anchors:
  - disease_stage_coarse == "Healthy" AND NOT dubious  =>  F0 (clean-healthy)
  - disease_stage_coarse == "Cirrhosis"                 =>  F4

These give the classifier real cross-cohort F0/F4 grounding without needing
documented Kleiner F-scores.

Training set (~146 donors):
  - 58 documented Andrews (F0=8, F1=8, F2=12, F3=11, F4=19)
  - 69 clean-healthy  -> F0
  - 28 Cirrhosis      -> F4

Excluded from training:
  - 35 dubious-healthy (ambiguous labels)
  - 101 Steatohepatitis (prediction targets; biologically F1-F3)
  - 36 Steatosis (held out for secondary biological check)

Outputs
-------
  - donor_fstage_augmented_predicted.tsv  (269 donors x predictions)
  - augmented_method_metrics.tsv          (head-to-head numbers)
  - figS_augmented_fstage_comparison.pdf  (4-panel)
  - appends rows to fstage_method_comparison.tsv
"""

from __future__ import annotations
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

from sklearn.metrics import cohen_kappa_score, confusion_matrix
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

try:
    import mord
    HAVE_MORD = True
except ImportError:
    HAVE_MORD = False

# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------
PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))

STAGE_DIR = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
FSTAGE_DOC    = STAGE_DIR / "donor_fstage_documented.tsv"
DUBIOUS       = STAGE_DIR / "dubious_healthy_donors.tsv"
DONOR_META    = STAGE_DIR / "donor_metadata_extended.tsv"
HEP_PB        = STAGE_DIR / "donor_hep_scvi_mean.tsv"
CELLFRAC      = STAGE_DIR / "cellfrac_by_stage.tsv"
METHOD_TBL    = STAGE_DIR / "fstage_method_comparison.tsv"
HEP_H5AD      = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_atlas_annotated.h5ad"

OUT_PRED      = STAGE_DIR / "donor_fstage_augmented_predicted.tsv"
OUT_METRICS   = STAGE_DIR / "augmented_method_metrics.tsv"
OUT_PDF       = PROJECT_ROOT / "figures/supplementary/stage_ccc/figS_augmented_fstage_comparison.pdf"

LATENT_COLS = [f"z{i}" for i in range(20)]


# ----------------------------------------------------------------------
# Color palette (publication-grade; matches masld_colors / control-gray rule)
# ----------------------------------------------------------------------
MASLD_PAL = {
    "Healthy":         "#9E9E9E",  # control gray (per CLAUDE.md rule)
    "Steatosis":       "#FFB454",
    "Steatohepatitis": "#E66F51",
    "Cirrhosis":       "#7A0E0E",
    "F0":              "#9E9E9E",
    "F1":              "#FFD16B",
    "F2":              "#FFB454",
    "F3":              "#E66F51",
    "F4":              "#7A0E0E",
    "vanilla":         "#5B7FB1",
    "augmented_ord":   "#2E8B57",
    "augmented_knn":   "#3D9970",
}


# ======================================================================
# Data assembly
# ======================================================================
def ensure_hep_pb() -> pd.DataFrame:
    """Load donor-level hepatocyte scVI mean from cache; recompute if missing."""
    if HEP_PB.exists():
        pb = pd.read_csv(HEP_PB, sep="\t")
        if all(c in pb.columns for c in LATENT_COLS):
            print(f"[load] hep pseudobulk cache: {len(pb)} donors")
            return pb
        print(f"[warn] hep pseudobulk cache missing latent columns; recomputing")

    print(f"[recompute] loading {HEP_H5AD.name} backed mode")
    a = ad.read_h5ad(HEP_H5AD, backed="r")
    if "X_scVI" not in a.obsm:
        sys.exit(f"FATAL: X_scVI missing from {HEP_H5AD}")
    Z = np.asarray(a.obsm["X_scVI"])
    samples = a.obs["sample"].astype(str).values
    datasets = a.obs["dataset"].astype(str).values

    df = pd.DataFrame(Z, columns=LATENT_COLS)
    df["sample"] = samples
    df["dataset"] = datasets
    pb = df.groupby(["sample", "dataset"], observed=True)[LATENT_COLS].mean().reset_index()
    n_hep = df.groupby(["sample", "dataset"], observed=True).size().rename(
        "n_hepatocytes_for_latent"
    ).reset_index()
    pb = pb.merge(n_hep, on=["sample", "dataset"])
    a.file.close()

    pb.to_csv(HEP_PB, sep="\t", index=False)
    print(f"[recompute] wrote cache {HEP_PB} ({len(pb)} donors)")
    return pb


def load_donors() -> pd.DataFrame:
    meta = pd.read_csv(DONOR_META, sep="\t")
    doc  = pd.read_csv(FSTAGE_DOC, sep="\t")
    doc["F_stage_documented"] = pd.to_numeric(doc["F_stage_documented"], errors="coerce")
    pb   = ensure_hep_pb()
    dub  = pd.read_csv(DUBIOUS, sep="\t")

    # Merge: meta (canonical 269) + latent + documented + dubious flag
    df = meta[["sample", "dataset", "disease_stage_coarse",
               "disease_stage_numeric"]].copy()
    df = df.merge(pb[["sample", "dataset"] + LATENT_COLS
                     + ["n_hepatocytes_for_latent"]],
                  on=["sample", "dataset"], how="left")
    df = df.merge(doc[["sample", "dataset", "F_stage_documented"]],
                  on=["sample", "dataset"], how="left")
    df["dubious"] = df["sample"].isin(dub["sample"].values)
    return df


def assemble_training(df: pd.DataFrame) -> pd.DataFrame:
    """Return training rows tagged with origin and y_train."""
    train = df.copy()
    train["y_train"] = np.nan
    train["origin"]  = "none"

    # 1. Documented Andrews
    mask_doc = train["F_stage_documented"].notna()
    train.loc[mask_doc, "y_train"] = train.loc[mask_doc, "F_stage_documented"]
    train.loc[mask_doc, "origin"]  = "documented"

    # 2. Clean-healthy => F0 (only if not already documented)
    mask_chealthy = ((train["disease_stage_coarse"] == "Healthy") &
                     (~train["dubious"]) &
                     (~mask_doc))
    train.loc[mask_chealthy, "y_train"] = 0
    train.loc[mask_chealthy, "origin"]  = "clean_healthy_anchor"

    # 3. Cirrhosis => F4 (only if not already documented)
    mask_cirr = ((train["disease_stage_coarse"] == "Cirrhosis") &
                 (~mask_doc))
    train.loc[mask_cirr, "y_train"] = 4
    train.loc[mask_cirr, "origin"]  = "cirrhosis_anchor"

    train_set = train[train["y_train"].notna() &
                      train[LATENT_COLS].notna().all(axis=1)].copy()
    train_set["y_train"] = train_set["y_train"].astype(int)
    return train_set


# ======================================================================
# Classifiers
# ======================================================================
def make_ordinal():
    if HAVE_MORD:
        return Pipeline([("scale", StandardScaler()),
                         ("clf",   mord.LogisticIT(alpha=1.0))])
    return Pipeline([("scale", StandardScaler()),
                     ("clf",   LogisticRegression(solver="lbfgs", C=1.0,
                                                   max_iter=2000,
                                                   class_weight="balanced"))])


def make_knn5():
    return Pipeline([("scale", StandardScaler()),
                     ("clf",   KNeighborsClassifier(n_neighbors=5,
                                                     weights="distance"))])


def predict_argmax(clf, X):
    if HAVE_MORD and isinstance(clf.named_steps["clf"], mord.LogisticIT):
        return clf.predict(X).astype(int)
    proba = clf.predict_proba(X)
    return clf.named_steps["clf"].classes_[proba.argmax(axis=1)].astype(int)


def predict_proba_full(clf, X):
    proba = clf.predict_proba(X)
    classes = clf.named_steps["clf"].classes_
    full = np.zeros((proba.shape[0], 5))
    for j, c in enumerate(classes):
        full[:, int(c)] = proba[:, j]
    rowsum = full.sum(axis=1, keepdims=True)
    rowsum[rowsum == 0] = 1.0
    return full / rowsum


# ======================================================================
# Cross-validation routines
# ======================================================================
def loocv(X, y, mk_clf, eval_mask=None):
    """LOOCV: hold one donor out at a time. eval_mask selects which
    held-out donors to score (default = all)."""
    n = len(y)
    y_pred = np.full(n, -1, dtype=int)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        clf = mk_clf()
        clf.fit(X[mask], y[mask])
        y_pred[i] = predict_argmax(clf, X[i:i+1])[0]
    if eval_mask is None:
        eval_mask = np.ones(n, dtype=bool)
    qwk = cohen_kappa_score(y[eval_mask], y_pred[eval_mask],
                            weights="quadratic", labels=list(range(5)))
    return float(qwk), y_pred


# ======================================================================
# Main pipeline
# ======================================================================
def run_method(name, mk_clf, df_all, train_set):
    """Run LOOCV variants + full predict for one classifier. Return dict."""
    print(f"\n========== METHOD: {name} ==========")
    X = train_set[LATENT_COLS].to_numpy()
    y = train_set["y_train"].to_numpy().astype(int)
    origin = train_set["origin"].to_numpy()

    # A. Full LOOCV across all 146 anchor+documented donors
    qwk_full, ypred_full = loocv(X, y, mk_clf)
    print(f"[A] Full-LOOCV (n={len(y)}): QWK={qwk_full:.4f}")
    cm_full = confusion_matrix(y, ypred_full, labels=list(range(5)))
    print("[A] CM rows=true F0..F4, cols=pred F0..F4")
    print(cm_full)

    # B. Andrews-only LOOCV (Andrews held out one at a time; anchors stay in
    #    training). This is apples-to-apples vs vanilla 0.742.
    andrews_idx = np.where(origin == "documented")[0]
    n = len(y)
    ypred_andrews = np.full(n, -1, dtype=int)
    for k, i in enumerate(andrews_idx):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        clf = mk_clf()
        clf.fit(X[mask], y[mask])
        ypred_andrews[i] = predict_argmax(clf, X[i:i+1])[0]
    qwk_andrews = cohen_kappa_score(
        y[andrews_idx], ypred_andrews[andrews_idx],
        weights="quadratic", labels=list(range(5)))
    print(f"[B] Andrews-LOOCV (n={len(andrews_idx)}): "
          f"QWK={qwk_andrews:.4f}  (vanilla 343b = 0.7422)")
    cm_and = confusion_matrix(y[andrews_idx], ypred_andrews[andrews_idx],
                              labels=list(range(5)))
    print("[B] Andrews CM:")
    print(cm_and)

    # F. Healthy->F4 hallucinations: clean-healthy donor LOOCV
    chealthy_idx = np.where(origin == "clean_healthy_anchor")[0]
    ypred_ch = np.full(n, -1, dtype=int)
    for i in chealthy_idx:
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        clf = mk_clf()
        clf.fit(X[mask], y[mask])
        ypred_ch[i] = predict_argmax(clf, X[i:i+1])[0]
    n_F4_hallucinate_clean = int((ypred_ch[chealthy_idx] == 4).sum())
    n_chealthy = len(chealthy_idx)
    print(f"[F] Clean-healthy LOOCV: {n_F4_hallucinate_clean}/{n_chealthy} "
          f"predicted F4 (vanilla all-Healthy = 27)")

    # Cirrhosis LOOCV for symmetry (F0 hallucinations from F4 anchors)
    cirr_idx = np.where(origin == "cirrhosis_anchor")[0]
    ypred_cirr = np.full(n, -1, dtype=int)
    for i in cirr_idx:
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        clf = mk_clf()
        clf.fit(X[mask], y[mask])
        ypred_cirr[i] = predict_argmax(clf, X[i:i+1])[0]
    n_F0_hallucinate_cirr = int((ypred_cirr[cirr_idx] == 0).sum())

    # Refit on full 146, predict everyone in df_all who has latents
    clf = mk_clf()
    clf.fit(X, y)
    pred_mask = df_all[LATENT_COLS].notna().all(axis=1)
    X_all = df_all.loc[pred_mask, LATENT_COLS].to_numpy()
    pred_all = predict_argmax(clf, X_all)
    proba_all = predict_proba_full(clf, X_all)
    pred_col = pd.Series(np.nan, index=df_all.index)
    pred_col.loc[pred_mask] = pred_all.astype(float)

    proba_cols = {}
    for k in range(5):
        s = pd.Series(np.nan, index=df_all.index)
        s.loc[pred_mask] = proba_all[:, k]
        proba_cols[f"P_F{k}_{name}"] = s

    # C. SH prediction breakdown
    sh_mask = df_all["disease_stage_coarse"] == "Steatohepatitis"
    sh_pred = pred_col[sh_mask & pred_mask].astype(int)
    sh_breakdown = sh_pred.value_counts(normalize=True).reindex(range(5),
                                                                 fill_value=0)
    sh_counts = sh_pred.value_counts().reindex(range(5), fill_value=0)
    print(f"[C] SH predictions (n={len(sh_pred)}):")
    for k in range(5):
        print(f"    F{k}: {int(sh_counts[k]):3d}  "
              f"({sh_breakdown[k]*100:5.1f}%)")
    print("    vanilla scVI cell-fraction SH baseline: "
          "F0=6.9%, F1=14.0%, F2=18.7%, F3=31.4%, F4=29.1%")

    # D. Steatosis breakdown
    st_mask = df_all["disease_stage_coarse"] == "Steatosis"
    st_pred = pred_col[st_mask & pred_mask].astype(int)
    st_breakdown = st_pred.value_counts(normalize=True).reindex(range(5),
                                                                 fill_value=0)
    st_counts = st_pred.value_counts().reindex(range(5), fill_value=0)
    print(f"[D] Steatosis predictions (n={len(st_pred)}):")
    for k in range(5):
        print(f"    F{k}: {int(st_counts[k]):3d}  "
              f"({st_breakdown[k]*100:5.1f}%)")

    # E. External Spearman vs disease_stage_numeric
    valid = df_all.assign(F_pred=pred_col).dropna(subset=[
        "F_pred", "disease_stage_numeric"])
    rho, p = stats.spearmanr(valid["F_pred"], valid["disease_stage_numeric"])
    print(f"[E] External Spearman rho vs disease_stage_numeric: "
          f"rho={rho:.4f}, p={p:.3e}, n={len(valid)}  "
          f"(vanilla 0.391)")

    # Per-dataset breakdown
    perds = []
    for ds, g in valid.groupby("dataset"):
        if len(g) >= 3:
            r, pp = stats.spearmanr(g["F_pred"], g["disease_stage_numeric"])
            perds.append((ds, len(g), r, pp))
    perds_df = pd.DataFrame(perds, columns=["dataset", "n", "rho", "p"])
    print(f"[E] Per-dataset rho:")
    print(perds_df.to_string(index=False))

    return dict(
        name=name,
        n_train=len(y),
        qwk_full_loocv=qwk_full,
        qwk_andrews_loocv=qwk_andrews,
        external_rho=float(rho),
        external_rho_p=float(p),
        n_external=len(valid),
        n_healthy_to_F4=n_F4_hallucinate_clean,
        n_chealthy=n_chealthy,
        n_F0_hallucinate_cirr=n_F0_hallucinate_cirr,
        n_cirr=len(cirr_idx),
        sh_counts=sh_counts.tolist(),
        sh_frac=sh_breakdown.round(4).tolist(),
        n_sh=len(sh_pred),
        st_counts=st_counts.tolist(),
        st_frac=st_breakdown.round(4).tolist(),
        n_st=len(st_pred),
        pred_col=pred_col,
        proba_cols=proba_cols,
        cm_andrews=cm_and,
        perds=perds_df,
    )


# ======================================================================
# Figure
# ======================================================================
def make_figure(results: list[dict], df_all: pd.DataFrame):
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    SH_BASELINE = [0.069, 0.140, 0.187, 0.314, 0.291]  # vanilla scVI cell-frac

    with PdfPages(OUT_PDF) as pdf:
        fig, axes = plt.subplots(2, 2, figsize=(13, 11))

        # Panel A: QWK & rho across methods
        ax = axes[0, 0]
        labels = ["vanilla\nscVI ordinal\n(343b)"] + \
                 [r["name"].replace("_", "\n") for r in results]
        x = np.arange(len(labels))
        w = 0.27
        andrews = [0.7422] + [r["qwk_andrews_loocv"] for r in results]
        full    = [np.nan]  + [r["qwk_full_loocv"]    for r in results]
        rho_ext = [0.391]   + [r["external_rho"]      for r in results]
        ax.bar(x - w, andrews, w, color=MASLD_PAL["vanilla"],
               label="Andrews-LOOCV QWK")
        ax.bar(x,     full,    w, color=MASLD_PAL["augmented_ord"],
               label="Full-LOOCV QWK")
        ax.bar(x + w, rho_ext, w, color=MASLD_PAL["augmented_knn"],
               label=r"External $\rho$ vs disease_stage")
        ax.axhline(0.7422, ls="--", lw=0.7, color="black", alpha=0.5)
        ax.axhline(0.391,  ls=":",  lw=0.7, color="black", alpha=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel("Score")
        ax.set_ylim(0, 1.0)
        ax.set_title("(a) Head-to-head: QWK and external Spearman")
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(axis="y", alpha=0.25)

        # Panel B: SH F-stage distribution (predicted) per method + baseline
        ax = axes[0, 1]
        stages = ["F0", "F1", "F2", "F3", "F4"]
        x = np.arange(5)
        w = 0.27
        ax.bar(x - w, SH_BASELINE, w, color=MASLD_PAL["vanilla"],
               label="vanilla scVI cell-frac")
        ax.bar(x,     results[0]["sh_frac"], w,
               color=MASLD_PAL["augmented_ord"],
               label=results[0]["name"])
        if len(results) > 1:
            ax.bar(x + w, results[1]["sh_frac"], w,
                   color=MASLD_PAL["augmented_knn"],
                   label=results[1]["name"])
        ax.set_xticks(x)
        ax.set_xticklabels(stages)
        ax.set_ylabel("Fraction of SH donors")
        ax.set_title(f"(b) Steatohepatitis donor F-stage "
                     f"(n={results[0]['n_sh']})")
        ax.legend(fontsize=8)
        ax.grid(axis="y", alpha=0.25)

        # Panel C: Steatosis F-stage distribution
        ax = axes[1, 0]
        ax.bar(x - w/2, results[0]["st_frac"], w,
               color=MASLD_PAL["augmented_ord"],
               label=results[0]["name"])
        if len(results) > 1:
            ax.bar(x + w/2, results[1]["st_frac"], w,
                   color=MASLD_PAL["augmented_knn"],
                   label=results[1]["name"])
        ax.set_xticks(x)
        ax.set_xticklabels(stages)
        ax.set_ylabel("Fraction of Steatosis donors")
        ax.set_title(f"(c) Steatosis donor F-stage held-out "
                     f"(n={results[0]['n_st']})")
        ax.legend(fontsize=8)
        ax.grid(axis="y", alpha=0.25)

        # Panel D: Healthy->F4 hallucinations
        ax = axes[1, 1]
        labels_d = ["vanilla\nscVI ordinal\n(343b)"] + \
                   [r["name"].replace("_", "\n") for r in results]
        vals = [27] + [r["n_healthy_to_F4"] for r in results]
        denoms = [69] + [r["n_chealthy"] for r in results]
        # vanilla scored ALL Healthy (104) without dubious filtering; show both
        bars = ax.bar(np.arange(len(labels_d)), vals,
                      color=[MASLD_PAL["vanilla"], MASLD_PAL["augmented_ord"],
                             MASLD_PAL["augmented_knn"]][:len(labels_d)])
        for i, (v, n) in enumerate(zip(vals, denoms)):
            ax.text(i, v + 0.5, f"{v}/{n}", ha="center", fontsize=9)
        ax.set_xticks(np.arange(len(labels_d)))
        ax.set_xticklabels(labels_d, fontsize=8)
        ax.set_ylabel("Healthy donors predicted F4")
        ax.set_title("(d) Hallucination check: clean-healthy LOOCV")
        ax.grid(axis="y", alpha=0.25)
        ax.set_ylim(0, max(vals) * 1.3 + 1)

        fig.suptitle(
            "Augmented F-stage classifier: anchor F0/F4 across cohorts",
            fontsize=12, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    print(f"[fig] wrote {OUT_PDF}")


# ======================================================================
# Outputs
# ======================================================================
def write_outputs(results, df_all, train_set):
    # Predictions table
    out = df_all[["sample", "dataset", "disease_stage_coarse",
                  "disease_stage_numeric"]].copy()
    out["F_stage_documented"] = df_all["F_stage_documented"]
    out["origin"] = df_all["sample"].map(
        train_set.set_index(["sample"])["origin"].to_dict()).fillna(
        "prediction_only")
    for r in results:
        out[f"F_stage_pred_{r['name']}"] = r["pred_col"]
        for col, vals in r["proba_cols"].items():
            out[col] = vals.round(6)
    out.to_csv(OUT_PRED, sep="\t", index=False)
    print(f"[output] {OUT_PRED}")

    # Metrics table
    rows = []
    for r in results:
        rows.append({
            "method": r["name"],
            "n_train": r["n_train"],
            "qwk_full_loocv": round(r["qwk_full_loocv"], 4),
            "qwk_andrews_loocv": round(r["qwk_andrews_loocv"], 4),
            "external_rho": round(r["external_rho"], 4),
            "external_rho_p": "{:.3e}".format(r["external_rho_p"]),
            "n_external": r["n_external"],
            "n_healthy_to_F4": r["n_healthy_to_F4"],
            "n_chealthy_loo": r["n_chealthy"],
            "n_cirrhosis_to_F0": r["n_F0_hallucinate_cirr"],
            "n_cirrhosis_loo": r["n_cirr"],
            "sh_pred_F0_frac": r["sh_frac"][0],
            "sh_pred_F1_frac": r["sh_frac"][1],
            "sh_pred_F2_frac": r["sh_frac"][2],
            "sh_pred_F3_frac": r["sh_frac"][3],
            "sh_pred_F4_frac": r["sh_frac"][4],
            "n_sh_pred":       r["n_sh"],
            "steat_pred_F0_frac": r["st_frac"][0],
            "steat_pred_F1_frac": r["st_frac"][1],
            "steat_pred_F2_frac": r["st_frac"][2],
            "steat_pred_F3_frac": r["st_frac"][3],
            "steat_pred_F4_frac": r["st_frac"][4],
            "n_steat_pred":       r["n_st"],
        })
    # Add vanilla reference row
    rows.append({
        "method": "vanilla_scvi_ordinal (343b ref)",
        "n_train": 58,
        "qwk_full_loocv": np.nan,
        "qwk_andrews_loocv": 0.7422,
        "external_rho": 0.391,
        "external_rho_p": "6.14e-11",
        "n_external": 260,
        "n_healthy_to_F4": 27,
        "n_chealthy_loo": 104,
        "n_cirrhosis_to_F0": 0,
        "n_cirrhosis_loo": 28,
        "sh_pred_F0_frac": 0.069,
        "sh_pred_F1_frac": 0.140,
        "sh_pred_F2_frac": 0.187,
        "sh_pred_F3_frac": 0.314,
        "sh_pred_F4_frac": 0.291,
        "n_sh_pred":       101,
        "steat_pred_F0_frac": np.nan,
        "steat_pred_F1_frac": np.nan,
        "steat_pred_F2_frac": np.nan,
        "steat_pred_F3_frac": np.nan,
        "steat_pred_F4_frac": np.nan,
        "n_steat_pred":       36,
    })
    pd.DataFrame(rows).to_csv(OUT_METRICS, sep="\t", index=False)
    print(f"[output] {OUT_METRICS}")

    # Append to global method comparison
    cmp_rows = []
    for r in results:
        cmp_rows.append({
            "method": r["name"],
            "n_train": r["n_train"],
            "qwk_loocv": round(r["qwk_andrews_loocv"], 4),
            "accuracy_loocv": "",
            "off_by_one_loocv": "",
        })
    cmp_df_new = pd.DataFrame(cmp_rows)

    if METHOD_TBL.exists():
        cmp_existing = pd.read_csv(METHOD_TBL, sep="\t")
        # Drop any prior augmented rows so we are idempotent
        cmp_existing = cmp_existing[~cmp_existing["method"].isin(
            cmp_df_new["method"].tolist())]
        cmp_all = pd.concat([cmp_existing, cmp_df_new], ignore_index=True)
    else:
        cmp_all = cmp_df_new
    cmp_all.to_csv(METHOD_TBL, sep="\t", index=False)
    print(f"[output] appended to {METHOD_TBL}")


# ======================================================================
def main():
    print(f"HAVE_MORD={HAVE_MORD}")
    OUT_PRED.parent.mkdir(parents=True, exist_ok=True)

    df_all = load_donors()
    print(f"[load] donor roster: {len(df_all)} donors")
    n_with_latent = int(df_all[LATENT_COLS].notna().all(axis=1).sum())
    print(f"[load] {n_with_latent}/{len(df_all)} donors have scVI latents")

    train_set = assemble_training(df_all)
    print(f"\n[train] augmented training set: n={len(train_set)}")
    print("[train] composition:")
    print(train_set.groupby(["origin", "y_train"]).size().unstack(fill_value=0))

    # Sanity: expect ~58 documented + 69 clean-healthy + 28 cirrhosis = 155 max
    # (Some may drop due to missing latents.)

    results = []
    results.append(run_method("augmented_scvi_ordinal",
                              make_ordinal, df_all, train_set))
    results.append(run_method("augmented_scvi_knn5",
                              make_knn5, df_all, train_set))

    write_outputs(results, df_all, train_set)
    make_figure(results, df_all)

    # Final summary
    print("\n" + "=" * 70)
    print("HEAD-TO-HEAD SUMMARY")
    print("=" * 70)
    print(f"{'method':<32s}  {'Andrews-LOOCV':>13s}  "
          f"{'full-LOOCV':>10s}  {'ext rho':>8s}  {'H->F4':>6s}")
    print(f"{'vanilla_scvi_ordinal (343b)':<32s}  "
          f"{0.7422:>13.4f}  {'--':>10s}  {0.391:>8.4f}  "
          f"{'27/104':>6s}")
    for r in results:
        print(f"{r['name']:<32s}  "
              f"{r['qwk_andrews_loocv']:>13.4f}  "
              f"{r['qwk_full_loocv']:>10.4f}  "
              f"{r['external_rho']:>8.4f}  "
              f"{r['n_healthy_to_F4']}/{r['n_chealthy']:<5d}")
    print("=" * 70)


if __name__ == "__main__":
    main()
