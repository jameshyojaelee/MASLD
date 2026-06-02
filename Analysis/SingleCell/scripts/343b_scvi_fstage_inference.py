#!/usr/bin/env python
"""
343b_scvi_fstage_inference.py

Project F-stage labels learned from 58 Andrews donors (GSE202379) onto the
remaining 211 donors using donor-level pseudobulk of the scVI latent space.

Why: F-stage is currently an Andrews-only axis. scVI integrated the atlas;
its latent space encodes biology common across datasets. If F-stage biology
is in the latent space, an Andrews-trained classifier should transfer.

Inputs
------
- Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_atlas_annotated.h5ad
    Contains obsm['X_scVI'] (657,804 hepatocytes x 20 dims) — the scVI latent
    used by the hepatocyte subclustering pipeline (Phase 1 trained scVI months
    ago; we do NOT retrain).
- Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_documented.tsv
    58 Andrews donors with Kleiner-scored F-stage 0-4.
- Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv
    269-donor roster.

Outputs
-------
- Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_scvi_predicted.tsv
    Columns: sample, dataset, F_stage_predicted_argmax, P_F0..P_F4,
             classifier_qwk_loocv, n_hepatocytes

Pipeline:
1. Build donor-level latent pseudobulk by averaging X_scVI across each
   donor's hepatocytes.
2. Train ordinal classifier (mord.LogisticIT) on the 58 Andrews donors.
3. LOOCV within Andrews -> QWK. Halt if < 0.15.
4. Refit on all 58, predict posteriors for all 269 donors.
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
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
try:
    import mord  # ordinal logistic
    HAVE_MORD = True
except ImportError:
    HAVE_MORD = False

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))

HEP_H5AD   = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_atlas_annotated.h5ad"
FSTAGE_DOC = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_documented.tsv"
DONOR_META = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv"
OUT_DIR    = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
OUT_TSV    = OUT_DIR / "donor_fstage_scvi_predicted.tsv"

QWK_HALT_THRESHOLD = 0.15  # below this, latent space lacks F-stage signal


def build_donor_pseudobulk() -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    """Return donor metadata table + latent pseudobulk matrix [n_donors x n_dim]."""
    print(f"[load] reading {HEP_H5AD.name} (backed)")
    a = ad.read_h5ad(HEP_H5AD, backed="r")
    if "X_scVI" not in a.obsm:
        sys.exit(f"FATAL: X_scVI missing from {HEP_H5AD}; obsm keys = {list(a.obsm.keys())}. "
                 "Per task spec, do NOT retrain scVI; halt.")
    Z = np.asarray(a.obsm["X_scVI"])  # 657,804 x 20
    samples = a.obs["sample"].astype(str).values
    datasets = a.obs["dataset"].astype(str).values
    print(f"[load] hepatocyte atlas: n_cells={Z.shape[0]}, n_latent_dim={Z.shape[1]}")

    df = pd.DataFrame({"sample": samples, "dataset": datasets})
    # group-wise mean
    # use pandas groupby on a DataFrame that carries Z columns
    cols = [f"z{i}" for i in range(Z.shape[1])]
    grouped = pd.DataFrame(Z, columns=cols)
    grouped["sample"] = samples
    grouped["dataset"] = datasets
    pb = grouped.groupby(["sample", "dataset"], observed=True)[cols].mean().reset_index()
    n_hep = grouped.groupby(["sample", "dataset"], observed=True).size().rename("n_hepatocytes").reset_index()
    pb = pb.merge(n_hep, on=["sample", "dataset"])
    a.file.close()
    print(f"[pseudobulk] {len(pb)} donor x dataset rows")
    return pb, cols


def attach_fstage(pb: pd.DataFrame) -> pd.DataFrame:
    fs = pd.read_csv(FSTAGE_DOC, sep="\t")
    fs["F_stage_documented"] = pd.to_numeric(fs["F_stage_documented"], errors="coerce")
    pb = pb.merge(
        fs[["sample", "dataset", "F_stage_documented"]],
        on=["sample", "dataset"], how="left",
    )
    return pb


def loocv_ordinal(X: np.ndarray, y: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """Leave-one-donor-out CV. Returns (QWK, y_true, y_pred)."""
    n = len(y)
    y_pred = np.empty(n, dtype=int)
    classes_sorted = np.sort(np.unique(y))
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        clf = make_classifier()
        clf.fit(X[mask], y[mask])
        # predict label by argmax of class posterior; map to integer scale
        if HAVE_MORD:
            y_pred[i] = int(clf.predict(X[i:i+1])[0])
        else:
            proba = clf.predict_proba(X[i:i+1])[0]
            y_pred[i] = int(clf.classes_[np.argmax(proba)])
    qwk = cohen_kappa_score(y, y_pred, weights="quadratic",
                             labels=list(range(5)))
    return float(qwk), y, y_pred


def make_classifier():
    """Ordinal classifier: prefer mord, fall back to multinomial logistic."""
    if HAVE_MORD:
        return Pipeline([
            ("scale", StandardScaler()),
            ("clf",   mord.LogisticIT(alpha=1.0)),
        ])
    # sklearn >=1.5 dropped multi_class kwarg; lbfgs auto-handles multinomial
    return Pipeline([
        ("scale", StandardScaler()),
        ("clf",   LogisticRegression(
            solver="lbfgs", C=1.0, max_iter=2000,
            class_weight="balanced",
        )),
    ])


def predict_with_proba(clf, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (argmax_label, proba_F0..F4).

    Mord LogisticIT has predict_proba in scikit-learn-style API (returns P over
    ordered classes). For LogisticRegression, .classes_ may not include all of
    {0..4} if some are absent from training; we pad zeros.
    """
    if HAVE_MORD:
        proba = clf.predict_proba(X)  # shape (n, n_classes_present)
        classes = clf.named_steps["clf"].classes_
    else:
        proba = clf.predict_proba(X)
        classes = clf.named_steps["clf"].classes_

    # Pad to full F0-F4 columns
    full = np.zeros((proba.shape[0], 5), dtype=float)
    for j, c in enumerate(classes):
        full[:, int(c)] = proba[:, j]
    # renormalize (rows sum to 1 only over present classes; fill missing as 0)
    rowsum = full.sum(axis=1, keepdims=True)
    rowsum[rowsum == 0] = 1.0
    full = full / rowsum
    argmax = full.argmax(axis=1)
    return argmax, full


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    pb, latent_cols = build_donor_pseudobulk()
    pb = attach_fstage(pb)

    donor = pd.read_csv(DONOR_META, sep="\t")
    print(f"[donor_meta] {len(donor)} donors in roster")

    # Use the donor roster as the ground truth list; restrict pseudobulk to
    # samples present in the roster (drops a handful of hep-atlas samples that
    # got excluded from donor_metadata for QC).
    pb = pb.merge(donor[["sample", "dataset"]], on=["sample", "dataset"], how="inner")
    print(f"[align] {len(pb)} donors with hepatocytes + in roster")

    train_mask = pb["F_stage_documented"].notna()
    n_train = int(train_mask.sum())
    print(f"[train] {n_train} Andrews donors with documented F-stage")

    if n_train < 20:
        sys.exit(f"FATAL: only {n_train} labeled donors; cannot train ordinal classifier.")

    X_train = pb.loc[train_mask, latent_cols].to_numpy()
    y_train = pb.loc[train_mask, "F_stage_documented"].astype(int).to_numpy()
    print(f"[train] X shape: {X_train.shape}; class counts: "
          f"{dict(zip(*np.unique(y_train, return_counts=True)))}")

    # LOOCV
    print("[loocv] running leave-one-donor-out within Andrews...")
    qwk, y_true, y_pred = loocv_ordinal(X_train, y_train)
    print(f"[loocv] QWK = {qwk:.3f}")
    cm = confusion_matrix(y_true, y_pred, labels=list(range(5)))
    print("[loocv] confusion matrix (rows=true F0..F4, cols=pred F0..F4):")
    print(cm)

    if qwk < QWK_HALT_THRESHOLD:
        sys.exit(f"HALT: QWK {qwk:.3f} < {QWK_HALT_THRESHOLD}; latent space does not "
                 "carry F-stage signal. Refusing to write predictions.")

    # Refit on all 58 Andrews, predict for everyone
    clf = make_classifier()
    clf.fit(X_train, y_train)
    X_all = pb[latent_cols].to_numpy()
    pred_argmax, pred_proba = predict_with_proba(clf, X_all)

    out = pb[["sample", "dataset", "n_hepatocytes"]].copy()
    out["F_stage_predicted_argmax"] = pred_argmax.astype(int)
    for k in range(5):
        out[f"P_F{k}"] = pred_proba[:, k].round(6)
    out["classifier_qwk_loocv"] = round(qwk, 4)

    # Ensure all 269 donors appear (even ones without hepatocytes -> NaN preds)
    full = donor[["sample", "dataset"]].merge(out, on=["sample", "dataset"], how="left")
    n_pred = full["F_stage_predicted_argmax"].notna().sum()
    print(f"[predict] {n_pred} / {len(full)} donors received predictions "
          f"(remainder lack hepatocytes)")

    full.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"[output] {OUT_TSV}")

    # Coverage summary
    print("\n[summary] Predicted F-stage distribution per dataset:")
    cov = (full.dropna(subset=["F_stage_predicted_argmax"])
                .assign(F=lambda d: d["F_stage_predicted_argmax"].astype(int))
                .groupby(["dataset", "F"], observed=True).size().unstack(fill_value=0))
    print(cov.to_string())


if __name__ == "__main__":
    main()
