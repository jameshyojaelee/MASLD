#!/usr/bin/env python
"""
343o_heldout_cohort_replication.py

Held-out cohort replication of the augmented F-stage classifier (Script 343m).

Per `docs/archive/single_cell_pre_remediation_2026-05-22/scrna_pipeline_rigor_plan.md`
Phase 7.4 (historical context in docs/technical/SINGLE_CELL_HISTORY.md) — the strongest external
validation. Hold out ALL donors from one dataset at a time, refit the
augmented scVI ordinal classifier on the remaining 5 cohorts' anchors, and
predict F-stage on the held-out dataset. Pre-registered threshold per
`docs/manuscript/pre_registered_thresholds.md`: AUROC > 0.65 on held-out
cohort = "generalizes". Below that = demote headline status of the
augmented method.

Procedure
---------
For each of 6 evaluable datasets (GSE202379, GSE244832, Liver_Atlas,
GSE185477, GSE136103, GSE174748 — exclude GSE189600 with n=2):
  1. Hold out ALL donors from that dataset.
  2. Refit augmented ordinal LogisticRegression on remaining 5 cohorts'
     anchors: 58 Andrews documented + remaining clean-healthy F0 +
     remaining cirrhosis F4.
  3. Predict F-stage on held-out cohort donors with latent.
  4. Compute:
     - AUROC for binary Healthy-vs-Diseased (true label = NOT Healthy in
       disease_stage_coarse). Score = P(F>=2) = P_F2 + P_F3 + P_F4.
     - Spearman rho vs disease_stage_numeric.
     - Confusion matrix vs disease_stage_coarse.
     - Healthy->F4 hallucinations among held-out healthy donors.

Andrews held-out is a special case: removing GSE202379 also removes ALL
documented F1/F2/F3 anchors, leaving only F0 (clean-healthy) and F4
(cirrhosis) endpoints. This tests whether endpoint anchors alone recover
the F1/F2/F3 ordinal structure.

Outputs
-------
  heldout_cohort_fstage_replication.tsv  per (held_out_cohort, donor)
  heldout_cohort_metrics.tsv             per held_out_cohort: AUROC/rho/etc
  figS_heldout_cohort_replication.pdf    6 small panels + AUROC bar chart
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import (
    roc_auc_score,
    confusion_matrix,
    accuracy_score,
    cohen_kappa_score,
)
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------
PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))

STAGE_DIR = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
FSTAGE_DOC = STAGE_DIR / "donor_fstage_documented.tsv"
DUBIOUS = STAGE_DIR / "dubious_healthy_donors.tsv"
DONOR_META = STAGE_DIR / "donor_metadata_extended.tsv"
HEP_PB = STAGE_DIR / "donor_hep_scvi_mean.tsv"

OUT_PRED = STAGE_DIR / "heldout_cohort_fstage_replication.tsv"
OUT_METRICS = STAGE_DIR / "heldout_cohort_metrics.tsv"
OUT_PDF = PROJECT_ROOT / "figures/supplementary/stage_ccc/figS_heldout_cohort_replication.pdf"

LATENT_COLS = [f"z{i}" for i in range(20)]

# Evaluable datasets (GSE189600 excluded with n=2)
EVAL_DATASETS = [
    "GSE202379",
    "GSE244832",
    "Liver_Atlas",
    "GSE185477",
    "GSE136103",
    "GSE174748",
]

# Pre-registered threshold from docs/manuscript/pre_registered_thresholds.md
AUROC_GENERALIZES_THRESHOLD = 0.65

# Color palette — match Script 343m / publication_theme.R / control-gray rule
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
    "pass":            "#2E8B57",
    "fail":            "#B22222",
}

STAGE_ORDER = ["Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"]


# ======================================================================
# Data assembly (parallels Script 343m)
# ======================================================================
def load_donors() -> pd.DataFrame:
    meta = pd.read_csv(DONOR_META, sep="\t")
    doc = pd.read_csv(FSTAGE_DOC, sep="\t")
    doc["F_stage_documented"] = pd.to_numeric(
        doc["F_stage_documented"], errors="coerce"
    )
    pb = pd.read_csv(HEP_PB, sep="\t")
    dub = pd.read_csv(DUBIOUS, sep="\t")

    df = meta[["sample", "dataset", "disease_stage_coarse",
               "disease_stage_numeric"]].copy()
    df = df.merge(
        pb[["sample", "dataset"] + LATENT_COLS],
        on=["sample", "dataset"], how="left"
    )
    df = df.merge(
        doc[["sample", "dataset", "F_stage_documented"]],
        on=["sample", "dataset"], how="left"
    )
    df["dubious"] = df["sample"].isin(dub["sample"].values)
    return df


def assemble_anchors(df: pd.DataFrame) -> pd.DataFrame:
    """Tag every donor with y_train (F-stage anchor) and origin label.

    Mirrors Script 343m's assemble_training:
      - documented Andrews => use F_stage_documented
      - clean-healthy (not dubious, not documented) => F0
      - cirrhosis (not documented) => F4
    """
    out = df.copy()
    out["y_train"] = np.nan
    out["origin"] = "none"

    mask_doc = out["F_stage_documented"].notna()
    out.loc[mask_doc, "y_train"] = out.loc[mask_doc, "F_stage_documented"]
    out.loc[mask_doc, "origin"] = "documented"

    mask_ch = ((out["disease_stage_coarse"] == "Healthy") &
               (~out["dubious"]) & (~mask_doc))
    out.loc[mask_ch, "y_train"] = 0
    out.loc[mask_ch, "origin"] = "clean_healthy_anchor"

    mask_cirr = ((out["disease_stage_coarse"] == "Cirrhosis") & (~mask_doc))
    out.loc[mask_cirr, "y_train"] = 4
    out.loc[mask_cirr, "origin"] = "cirrhosis_anchor"

    return out


# ======================================================================
# Classifier
# ======================================================================
def make_ordinal():
    """Multinomial LR as ordinal proxy (mord unavailable in env).

    Note: sklearn >= 1.7 dropped the `multi_class` arg and defaults to
    multinomial for non-binary problems. We rely on that default.
    """
    return Pipeline([
        ("scale", StandardScaler()),
        ("clf", LogisticRegression(
            solver="lbfgs",
            C=1.0,
            max_iter=5000,
            class_weight="balanced",
        )),
    ])


def predict_proba_full(clf, X):
    """Return (n,5) prob matrix indexed by F0..F4."""
    proba = clf.predict_proba(X)
    classes = clf.named_steps["clf"].classes_
    full = np.zeros((proba.shape[0], 5))
    for j, c in enumerate(classes):
        full[:, int(c)] = proba[:, j]
    rowsum = full.sum(axis=1, keepdims=True)
    rowsum[rowsum == 0] = 1.0
    return full / rowsum


def predict_argmax_from_proba(proba_full):
    return proba_full.argmax(axis=1).astype(int)


# ======================================================================
# Per-cohort held-out evaluation
# ======================================================================
def evaluate_heldout(df_anchored: pd.DataFrame, holdout_ds: str) -> dict:
    """Train on all anchors NOT in holdout_ds, predict on holdout_ds donors."""
    print(f"\n========== HELD OUT: {holdout_ds} ==========")

    # Training set: anchors from datasets other than holdout
    train_mask = (
        df_anchored["y_train"].notna() &
        df_anchored[LATENT_COLS].notna().all(axis=1) &
        (df_anchored["dataset"] != holdout_ds)
    )
    train = df_anchored.loc[train_mask].copy()
    train["y_train"] = train["y_train"].astype(int)

    # Held-out set: everyone in holdout_ds with a latent
    holdout_mask = (
        (df_anchored["dataset"] == holdout_ds) &
        df_anchored[LATENT_COLS].notna().all(axis=1)
    )
    holdout = df_anchored.loc[holdout_mask].copy()

    print(f"[train] n={len(train)} from {train['dataset'].nunique()} cohorts")
    print(f"[train] origin breakdown:")
    print(train.groupby(["origin", "y_train"]).size().unstack(fill_value=0))
    print(f"[holdout] n={len(holdout)} from {holdout_ds}")

    if len(train) < 10 or len(holdout) == 0:
        print(f"[skip] insufficient data")
        return None

    # Andrews held-out special case: only F0/F4 in training (no F1/F2/F3)
    classes_train = sorted(train["y_train"].unique().tolist())
    print(f"[train] classes present: {classes_train}")
    is_endpoint_only = set(classes_train).issubset({0, 4})
    if is_endpoint_only:
        print(f"[NOTE] endpoint-only training (F0/F4); F1/F2/F3 not represented")

    # Fit
    X_tr = train[LATENT_COLS].to_numpy()
    y_tr = train["y_train"].to_numpy().astype(int)
    X_ho = holdout[LATENT_COLS].to_numpy()

    clf = make_ordinal()
    clf.fit(X_tr, y_tr)
    proba = predict_proba_full(clf, X_ho)
    pred = predict_argmax_from_proba(proba)

    holdout = holdout.assign(
        F_stage_pred=pred,
        P_F0=proba[:, 0], P_F1=proba[:, 1], P_F2=proba[:, 2],
        P_F3=proba[:, 3], P_F4=proba[:, 4],
    )

    # ---- Binary AUROC: Healthy vs Diseased ----
    y_binary = (holdout["disease_stage_coarse"] != "Healthy").astype(int)
    score_diseased = holdout["P_F2"] + holdout["P_F3"] + holdout["P_F4"]
    auroc = np.nan
    if y_binary.nunique() == 2:
        auroc = roc_auc_score(y_binary, score_diseased)
        print(f"[AUROC] Healthy-vs-Diseased: {auroc:.4f}  "
              f"(n_pos={int(y_binary.sum())}, n_neg={int((1-y_binary).sum())})")
    else:
        print(f"[AUROC] N/A (only one class present in holdout)")

    # ---- Spearman rho vs disease_stage_numeric ----
    rho, rho_p, n_rho = np.nan, np.nan, 0
    valid_rho = holdout.dropna(subset=["disease_stage_numeric"])
    if len(valid_rho) >= 3 and valid_rho["disease_stage_numeric"].nunique() >= 2:
        rho, rho_p = stats.spearmanr(
            valid_rho["F_stage_pred"], valid_rho["disease_stage_numeric"]
        )
        n_rho = len(valid_rho)
        print(f"[rho]   Spearman vs disease_stage_numeric: "
              f"rho={rho:.4f}, p={rho_p:.3e}, n={n_rho}")
    else:
        print(f"[rho]   N/A (insufficient variation)")

    # ---- Accuracy vs disease_stage_coarse (map pred F-stage to coarse) ----
    # Mapping: F0 -> Healthy; F1 -> Steatosis; F2,F3 -> SH; F4 -> Cirrhosis
    coarse_map = {0: "Healthy", 1: "Steatosis", 2: "Steatohepatitis",
                  3: "Steatohepatitis", 4: "Cirrhosis"}
    holdout["pred_coarse"] = holdout["F_stage_pred"].map(coarse_map)
    valid_acc = holdout.dropna(subset=["disease_stage_coarse"])
    acc = accuracy_score(valid_acc["disease_stage_coarse"],
                         valid_acc["pred_coarse"]) if len(valid_acc) else np.nan
    print(f"[acc]   coarse-class accuracy: {acc:.4f}  (n={len(valid_acc)})")

    # ---- Confusion matrix (rows = true coarse, cols = pred coarse) ----
    cm = pd.crosstab(
        valid_acc["disease_stage_coarse"], valid_acc["pred_coarse"],
        rownames=["true"], colnames=["pred"], dropna=False
    ).reindex(index=STAGE_ORDER, columns=STAGE_ORDER, fill_value=0)
    print("[CM]    rows=true coarse, cols=pred coarse:")
    print(cm)

    # ---- Healthy -> F4 hallucinations among held-out healthy donors ----
    h_mask = holdout["disease_stage_coarse"] == "Healthy"
    n_h = int(h_mask.sum())
    n_h_to_F4 = int(((holdout.loc[h_mask, "F_stage_pred"] == 4)).sum())
    n_h_to_diseased = int(
        (holdout.loc[h_mask, "F_stage_pred"] >= 2).sum()
    )
    print(f"[hallu] Healthy donors: {n_h_to_F4}/{n_h} predicted F4; "
          f"{n_h_to_diseased}/{n_h} predicted F>=2 (diseased)")

    # ---- Documented overlap: QWK vs documented F-stage if available ----
    # Important for Andrews held-out (the endpoint-only case): the only
    # apples-to-apples score with within-Andrews LOOCV QWK = 0.720.
    doc_mask = holdout["F_stage_documented"].notna()
    n_doc = int(doc_mask.sum())
    qwk_documented = np.nan
    if n_doc >= 3:
        qwk_documented = cohen_kappa_score(
            holdout.loc[doc_mask, "F_stage_documented"].astype(int),
            holdout.loc[doc_mask, "F_stage_pred"].astype(int),
            weights="quadratic", labels=list(range(5)),
        )
        print(f"[QWK]   vs documented F-stage: QWK={qwk_documented:.4f} "
              f"(n={n_doc})  -- compare to within-Andrews LOOCV 0.720")

    # ---- Healthy specificity: for Healthy-only holdouts, this is the
    #      only meaningful generalization metric ----
    healthy_specificity = np.nan
    if n_h > 0:
        # Specificity = fraction of true Healthy correctly classified as F<2
        n_h_correct = int((holdout.loc[h_mask, "F_stage_pred"] < 2).sum())
        healthy_specificity = n_h_correct / n_h
        print(f"[spec]  Healthy specificity (predicted F<2): "
              f"{healthy_specificity:.4f}  ({n_h_correct}/{n_h})")

    # ---- Pass/fail rule with edge-case handling ----
    # AUROC undefined when holdout has only one binary class.
    # In that case we fall back to:
    #   - QWK vs documented F-stage (if available; threshold 0.5 = halfway
    #     to within-Andrews 0.720)
    #   - Healthy specificity (if all-healthy holdout; threshold 0.5 = better
    #     than chance binary classification)
    if not np.isnan(auroc):
        passes_threshold = bool(auroc > AUROC_GENERALIZES_THRESHOLD)
        pass_metric = "AUROC"
    elif not np.isnan(qwk_documented):
        passes_threshold = bool(qwk_documented > 0.5)
        pass_metric = "QWK_doc"
    elif not np.isnan(healthy_specificity):
        passes_threshold = bool(healthy_specificity > 0.5)
        pass_metric = "spec_healthy"
    else:
        passes_threshold = False
        pass_metric = "undefined"

    return dict(
        holdout_dataset=holdout_ds,
        n_train=len(train),
        n_train_cohorts=int(train["dataset"].nunique()),
        endpoint_only_training=is_endpoint_only,
        n_holdout=len(holdout),
        n_holdout_healthy=n_h,
        auroc_hvd=float(auroc) if not np.isnan(auroc) else np.nan,
        spearman_rho=float(rho) if not np.isnan(rho) else np.nan,
        spearman_p=float(rho_p) if not np.isnan(rho_p) else np.nan,
        n_for_rho=n_rho,
        coarse_accuracy=float(acc) if not np.isnan(acc) else np.nan,
        n_for_accuracy=len(valid_acc),
        n_healthy_to_F4=n_h_to_F4,
        n_healthy_to_diseased=n_h_to_diseased,
        healthy_specificity=(
            float(healthy_specificity)
            if not np.isnan(healthy_specificity) else np.nan
        ),
        n_documented_in_holdout=n_doc,
        qwk_vs_documented=(
            float(qwk_documented)
            if not np.isnan(qwk_documented) else np.nan
        ),
        passes_threshold=passes_threshold,
        pass_metric=pass_metric,
        cm=cm,
        holdout_df=holdout,
    )


# ======================================================================
# Figure
# ======================================================================
def make_figure(results: list[dict]):
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    n = len(results)

    with PdfPages(OUT_PDF) as pdf:
        # ---- Page 1: per-cohort small panels (3 cols x 2 rows) ----
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()
        for i, r in enumerate(results):
            ax = axes[i]
            cm = r["cm"]
            im = ax.imshow(cm.values, cmap="Reds", aspect="auto")
            ax.set_xticks(range(len(STAGE_ORDER)))
            ax.set_xticklabels(STAGE_ORDER, rotation=35, ha="right", fontsize=8)
            ax.set_yticks(range(len(STAGE_ORDER)))
            ax.set_yticklabels(STAGE_ORDER, fontsize=8)
            ax.set_xlabel("Predicted coarse", fontsize=9)
            ax.set_ylabel("True coarse", fontsize=9)
            for ii in range(cm.shape[0]):
                for jj in range(cm.shape[1]):
                    val = cm.values[ii, jj]
                    if val > 0:
                        ax.text(jj, ii, str(int(val)), ha="center", va="center",
                                fontsize=8,
                                color="white" if val > cm.values.max() / 2
                                else "black")
            auroc_s = (f"{r['auroc_hvd']:.3f}"
                       if not np.isnan(r["auroc_hvd"]) else "NA")
            rho_s = (f"{r['spearman_rho']:.3f}"
                     if not np.isnan(r["spearman_rho"]) else "NA")
            qwk_s = (f"{r['qwk_vs_documented']:.3f}"
                     if not np.isnan(r["qwk_vs_documented"]) else "NA")
            pass_s = ("PASS" if r["passes_threshold"]
                      else ("FAIL" if r["pass_metric"] != "undefined"
                            else "N/A"))
            color = (MASLD_PAL["pass"] if r["passes_threshold"]
                     else MASLD_PAL["fail"])
            ep_tag = (" (endpoint-only train)" if r["endpoint_only_training"]
                      else "")
            metric_tag = f"[{r['pass_metric']}={pass_s}]"
            ax.set_title(
                f"({chr(97+i)}) {r['holdout_dataset']}{ep_tag}\n"
                f"AUROC={auroc_s}, rho={rho_s}, QWK_doc={qwk_s} "
                f"{metric_tag}\n"
                f"H->F4={r['n_healthy_to_F4']}/{r['n_holdout_healthy']}, "
                f"n={r['n_holdout']}",
                fontsize=8, color=color,
            )
        for j in range(len(results), len(axes)):
            axes[j].axis("off")
        fig.suptitle(
            "Held-out cohort F-stage replication "
            "(augmented scVI ordinal, leave-one-dataset-out)",
            fontsize=12, fontweight="bold",
        )
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # ---- Page 2: per-cohort AUROC bar chart with threshold line ----
        fig, ax = plt.subplots(figsize=(9, 5.5))
        names = [r["holdout_dataset"] for r in results]
        aurocs = [r["auroc_hvd"] if not np.isnan(r["auroc_hvd"]) else np.nan
                  for r in results]
        # Show AUROC where defined; for cohorts without AUROC use fallback metric
        plot_vals = []
        plot_labels = []
        plot_colors = []
        for r in results:
            if not np.isnan(r["auroc_hvd"]):
                plot_vals.append(r["auroc_hvd"])
                plot_labels.append(f"AUROC={r['auroc_hvd']:.3f}")
            elif not np.isnan(r["qwk_vs_documented"]):
                plot_vals.append(r["qwk_vs_documented"])
                plot_labels.append(f"QWK_doc={r['qwk_vs_documented']:.3f}")
            elif not np.isnan(r["healthy_specificity"]):
                plot_vals.append(r["healthy_specificity"])
                plot_labels.append(f"spec={r['healthy_specificity']:.3f}")
            else:
                plot_vals.append(0)
                plot_labels.append("NA")
            plot_colors.append(
                MASLD_PAL["pass"] if r["passes_threshold"]
                else MASLD_PAL["fail"]
            )
        bars = ax.bar(np.arange(len(names)), plot_vals, color=plot_colors,
                      edgecolor="black", linewidth=0.5)
        for i, b in enumerate(bars):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.012,
                    plot_labels[i], ha="center", fontsize=8)
        ax.axhline(AUROC_GENERALIZES_THRESHOLD, ls="--", color="black",
                   lw=1.0, alpha=0.65,
                   label=f"pre-registered AUROC threshold = "
                         f"{AUROC_GENERALIZES_THRESHOLD}")
        ax.axhline(0.5, ls=":", color="gray", lw=0.6,
                   label="random (0.5)")
        ax.set_xticks(np.arange(len(names)))
        ax.set_xticklabels(names, rotation=35, ha="right", fontsize=9)
        ax.set_ylabel("Held-out metric (AUROC, or fallback)")
        ax.set_ylim(0, 1.10)
        n_pass = sum(r["passes_threshold"] for r in results)
        ax.set_title(
            f"Per-cohort held-out replication\n"
            f"AUROC primary; fallback metrics for cohorts without "
            f"two-class holdout. {n_pass}/{len(results)} pass.",
            fontsize=10,
        )
        ax.legend(fontsize=8, loc="lower right")
        ax.grid(axis="y", alpha=0.25)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    print(f"[fig] wrote {OUT_PDF}")


# ======================================================================
# Outputs
# ======================================================================
def write_outputs(results: list[dict]):
    # Predictions (per donor x per held-out cohort)
    pred_rows = []
    for r in results:
        h = r["holdout_df"].copy()
        h["held_out_cohort"] = r["holdout_dataset"]
        h["endpoint_only_training"] = r["endpoint_only_training"]
        pred_rows.append(h[[
            "held_out_cohort", "sample", "dataset",
            "disease_stage_coarse", "disease_stage_numeric",
            "F_stage_documented", "F_stage_pred",
            "P_F0", "P_F1", "P_F2", "P_F3", "P_F4",
            "endpoint_only_training",
        ]])
    pred_df = pd.concat(pred_rows, ignore_index=True)
    pred_df.to_csv(OUT_PRED, sep="\t", index=False)
    print(f"[output] {OUT_PRED}  ({len(pred_df)} rows)")

    # Metrics
    metric_rows = []
    for r in results:
        metric_rows.append({
            "held_out_cohort": r["holdout_dataset"],
            "n_train": r["n_train"],
            "n_train_cohorts": r["n_train_cohorts"],
            "endpoint_only_training": r["endpoint_only_training"],
            "n_holdout": r["n_holdout"],
            "n_holdout_healthy": r["n_holdout_healthy"],
            "n_documented_in_holdout": r["n_documented_in_holdout"],
            "auroc_healthy_vs_diseased": (
                round(r["auroc_hvd"], 4)
                if not np.isnan(r["auroc_hvd"]) else np.nan
            ),
            "spearman_rho": (
                round(r["spearman_rho"], 4)
                if not np.isnan(r["spearman_rho"]) else np.nan
            ),
            "spearman_p": (
                "{:.3e}".format(r["spearman_p"])
                if not np.isnan(r["spearman_p"]) else "NA"
            ),
            "n_for_rho": r["n_for_rho"],
            "coarse_accuracy": (
                round(r["coarse_accuracy"], 4)
                if not np.isnan(r["coarse_accuracy"]) else np.nan
            ),
            "n_for_accuracy": r["n_for_accuracy"],
            "qwk_vs_documented": (
                round(r["qwk_vs_documented"], 4)
                if not np.isnan(r["qwk_vs_documented"]) else np.nan
            ),
            "healthy_specificity": (
                round(r["healthy_specificity"], 4)
                if not np.isnan(r["healthy_specificity"]) else np.nan
            ),
            "n_healthy_to_F4_hallucinations": r["n_healthy_to_F4"],
            "n_healthy_to_diseased": r["n_healthy_to_diseased"],
            "pass_metric": r["pass_metric"],
            "passes_threshold": r["passes_threshold"],
        })
    metric_df = pd.DataFrame(metric_rows)
    metric_df.to_csv(OUT_METRICS, sep="\t", index=False)
    print(f"[output] {OUT_METRICS}")


# ======================================================================
def main():
    OUT_PRED.parent.mkdir(parents=True, exist_ok=True)
    print(f"[info] PROJECT_ROOT = {PROJECT_ROOT}")
    print(f"[info] STAGE_DIR    = {STAGE_DIR}")
    print(f"[info] pre-registered AUROC threshold = "
          f"{AUROC_GENERALIZES_THRESHOLD}")

    df = load_donors()
    print(f"[load] donor roster: {len(df)} donors")
    n_with_latent = int(df[LATENT_COLS].notna().all(axis=1).sum())
    print(f"[load] {n_with_latent}/{len(df)} donors have scVI latents")

    df_anch = assemble_anchors(df)
    n_anch = int(df_anch["y_train"].notna().sum())
    print(f"[anchors] {n_anch} anchor donors total")
    print(df_anch.groupby(["origin", "y_train"]).size()
          .unstack(fill_value=0))
    print("[anchors] anchors per dataset:")
    print(df_anch[df_anch["y_train"].notna()].groupby(
        ["dataset", "origin"]).size().unstack(fill_value=0))

    # Run all 6 evaluable held-outs
    results = []
    for ds in EVAL_DATASETS:
        r = evaluate_heldout(df_anch, ds)
        if r is not None:
            results.append(r)

    write_outputs(results)
    make_figure(results)

    # Summary
    print("\n" + "=" * 90)
    print("HELD-OUT COHORT REPLICATION SUMMARY")
    print("=" * 90)
    print(f"{'cohort':<14s} {'n_train':>7s} {'n_hold':>6s} "
          f"{'AUROC':>6s} {'rho':>7s} {'QWK_doc':>7s} "
          f"{'H->F4':>7s} {'metric':>13s} {'PASS':>5s}")
    n_pass = 0
    n_eval = 0  # cohorts with a defined pass metric
    for r in results:
        passing = r["passes_threshold"]
        n_pass += int(passing)
        if r["pass_metric"] != "undefined":
            n_eval += 1
        au = (f"{r['auroc_hvd']:.3f}"
              if not np.isnan(r["auroc_hvd"]) else "   NA")
        rho = (f"{r['spearman_rho']:+.3f}"
               if not np.isnan(r["spearman_rho"]) else "   NA  ")
        qwk = (f"{r['qwk_vs_documented']:.3f}"
               if not np.isnan(r["qwk_vs_documented"]) else "   NA ")
        hh = f"{r['n_healthy_to_F4']}/{r['n_holdout_healthy']}"
        verdict = "PASS" if passing else (
            "FAIL" if r["pass_metric"] != "undefined" else " N/A"
        )
        print(f"{r['holdout_dataset']:<14s} "
              f"{r['n_train']:>7d} {r['n_holdout']:>6d} "
              f"{au:>6s} {rho:>7s} {qwk:>7s} {hh:>7s} "
              f"{r['pass_metric']:>13s} {verdict:>5s}")
    print("-" * 90)
    print(f"Cohorts passing (primary AUROC>0.65 or fallback): "
          f"{n_pass}/{len(results)} ({n_pass}/{n_eval} evaluable)")
    print(f"Reference: 343m within-Andrews LOOCV QWK = 0.720; "
          f"external rho = 0.639")
    if n_pass >= len(results) - 1:
        print("VERDICT: augmented method GENERALIZES to held-out cohorts "
              "(>=5/6 pass).")
    elif n_pass >= max(1, len(results) // 2):
        print("VERDICT: partial generalization -- keep with caveats; "
              "report per-cohort metrics in supplement.")
    else:
        print(f"VERDICT: DEMOTE -- pre-registered threshold violated in "
              f"{len(results) - n_pass}/{len(results)} cohorts.")
    print("=" * 90)


if __name__ == "__main__":
    main()
