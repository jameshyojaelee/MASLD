#!/usr/bin/env python
"""
343c_fstage_alt_inference_compare.py

Head-to-head comparison of alternative F-stage inference methods against the
scVI ordinal classifier (343b).

Methods (all trained on the 58 GSE202379/Andrews donors with documented F-stage):
  1. Cell-type proportion (9 lineage fractions) + ordinal logistic
  2. cNMF program scores (k=16 global) + ordinal logistic
  3. k-NN (k=5, median) in donor-mean hepatocyte scVI latent space
  4. Pseudotime axis (mac_pt + hep_pt + progressor_frac) + ordinal logistic

Outputs
-------
- fstage_method_comparison.tsv  -- per-method LOOCV QWK / accuracy / off-by-one
- donor_fstage_alt_predictions.tsv -- per-donor predictions for all 269
- fstage_method_confusion.txt -- per-method LOOCV confusion matrices

Constraints
-----------
- mord absent in rapids_singlecell -> uses sklearn LogisticRegression as
  ordinal proxy (same fallback the scVI agent uses).
- Does NOT retrain scVI. Recomputes donor-mean hepatocyte latent from
  hepatocyte_atlas_annotated.h5ad backed-read.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import cohen_kappa_score, confusion_matrix

try:
    import mord  # noqa: F401
    HAVE_MORD = True
except ImportError:
    HAVE_MORD = False

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
ST_DIR     = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
HEP_H5AD   = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_atlas_annotated.h5ad"
LINEAGE    = ST_DIR / "lineage_cell_counts.tsv"
FSTAGE_DOC = ST_DIR / "donor_fstage_documented.tsv"
EXT_META   = ST_DIR / "donor_metadata_extended.tsv"
DONOR_META = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv"
CNMF_WIDE  = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/integration/program_donor_scores_wide.tsv.gz"
SCVI_PRED  = ST_DIR / "donor_fstage_scvi_predicted.tsv"

OUT_COMP   = ST_DIR / "fstage_method_comparison.tsv"
OUT_PRED   = ST_DIR / "donor_fstage_alt_predictions.tsv"
OUT_CMTXT  = ST_DIR / "fstage_method_confusion.txt"
HEP_MEAN_CACHE = ST_DIR / "donor_hep_scvi_mean.tsv"

# ---------------------------------------------------------------------------
# Lineage mapping: collapse 16 granular cell types -> 9 lineages requested.
LINEAGE_MAP = {
    "Hep":    ["n_Hepatocytes"],
    "Mac":    ["n_Macrophages", "n_Mono+mono_derived_cells",
               "n_cDC1s", "n_cDC2s", "n_pDCs",
               "n_Basophils", "n_Neutrophils"],
    "Fib":    ["n_Fibroblasts"],
    "LSEC":   ["n_Endothelial_cells"],
    "Chol":   ["n_Cholangiocytes"],
    "T":      ["n_T_cells"],
    "B":      ["n_B_cells"],
    "NK":     ["n_Resident_NK", "n_Circulating_NK/NKT"],
    "Plasma": ["n_Plasma_cells"],
}


def make_ordinal_clf():
    if HAVE_MORD:
        import mord
        return Pipeline([("scale", StandardScaler()),
                         ("clf", mord.LogisticIT(alpha=1.0))])
    # sklearn fallback (same as scVI agent)
    return Pipeline([("scale", StandardScaler()),
                     ("clf", LogisticRegression(solver="lbfgs",
                                                C=1.0,
                                                max_iter=2000,
                                                class_weight="balanced"))])


def predict_with_full_proba(clf, X):
    """Return argmax label + (n, 5) proba matrix padded to F0..F4."""
    proba = clf.predict_proba(X)
    classes = clf.named_steps["clf"].classes_
    full = np.zeros((proba.shape[0], 5), dtype=float)
    for j, c in enumerate(classes):
        full[:, int(c)] = proba[:, j]
    rs = full.sum(axis=1, keepdims=True)
    rs[rs == 0] = 1.0
    full = full / rs
    return full.argmax(axis=1), full


def loocv_ordinal(X, y):
    n = len(y)
    y_pred = np.empty(n, dtype=int)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        clf = make_ordinal_clf()
        clf.fit(X[mask], y[mask])
        ar, _ = predict_with_full_proba(clf, X[i:i+1])
        y_pred[i] = int(ar[0])
    return y, y_pred


def loocv_knn_median(X, y, k=5):
    """k-NN regression by neighbor median; round to nearest integer F-stage."""
    n = len(y)
    y_pred = np.empty(n, dtype=int)
    # Standardize across train fold (refit per-fold for cleanliness;
    # X is donor-mean latent -> scaling is mild)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        Xt = X[mask]
        yt = y[mask]
        # use Euclidean on raw latent (it's already a learned embedding)
        d = np.linalg.norm(Xt - X[i], axis=1)
        idx = np.argpartition(d, kth=min(k, len(d)-1))[:k]
        pred = int(np.round(np.median(yt[idx])))
        pred = max(0, min(4, pred))
        y_pred[i] = pred
    return y, y_pred


def knn_predict_all(X_train, y_train, X_all, k=5):
    """Predict label + pseudo-posterior for all donors."""
    out_lab = np.empty(X_all.shape[0], dtype=int)
    out_proba = np.zeros((X_all.shape[0], 5), dtype=float)
    for i in range(X_all.shape[0]):
        d = np.linalg.norm(X_train - X_all[i], axis=1)
        idx = np.argpartition(d, kth=min(k, len(d)-1))[:k]
        neigh = y_train[idx]
        med = int(np.round(np.median(neigh)))
        med = max(0, min(4, med))
        out_lab[i] = med
        # empirical neighbor-frequency posterior
        for c in neigh:
            out_proba[i, int(c)] += 1
        out_proba[i] /= max(1, out_proba[i].sum())
    return out_lab, out_proba


def qwk_metrics(y_true, y_pred):
    qwk = cohen_kappa_score(y_true, y_pred, weights="quadratic",
                            labels=list(range(5)))
    acc = float(np.mean(y_true == y_pred))
    off1 = float(np.mean(np.abs(y_true - y_pred) <= 1))
    cm = confusion_matrix(y_true, y_pred, labels=list(range(5)))
    return float(qwk), acc, off1, cm


# ---------------------------------------------------------------------------
# Feature builders
# ---------------------------------------------------------------------------
def build_celltype_features(roster: pd.DataFrame) -> pd.DataFrame:
    """Return donor-level 9-lineage fraction matrix."""
    lc = pd.read_csv(LINEAGE, sep="\t")
    # collapse to 9 lineages
    out = pd.DataFrame({"sample": lc["sample"]})
    for lin, cols in LINEAGE_MAP.items():
        present = [c for c in cols if c in lc.columns]
        out[f"n_{lin}"] = lc[present].sum(axis=1) if present else 0
    totals = out[[c for c in out.columns if c.startswith("n_")]].sum(axis=1)
    for lin in LINEAGE_MAP:
        out[f"frac_{lin}"] = np.where(totals > 0, out[f"n_{lin}"] / totals, 0.0)
    out = out[["sample"] + [f"frac_{lin}" for lin in LINEAGE_MAP]]
    out = roster.merge(out, on="sample", how="left")
    return out


def build_cnmf_features(roster: pd.DataFrame) -> pd.DataFrame:
    cnmf = pd.read_csv(CNMF_WIDE, sep="\t")
    feat_cols = [c for c in cnmf.columns if c.startswith("cnmf_global_k16_")]
    out = roster.merge(cnmf[["sample"] + feat_cols], on="sample", how="left")
    return out


def build_pseudotime_features(roster: pd.DataFrame) -> pd.DataFrame:
    em = pd.read_csv(EXT_META, sep="\t",
                     usecols=["sample", "dataset",
                              "macrophage_pseudotime_mean",
                              "hepatocyte_pseudotime_mean",
                              "progressor_frac"])
    out = roster.merge(em, on=["sample", "dataset"], how="left")
    return out


def build_hep_scvi_mean(roster: pd.DataFrame) -> pd.DataFrame:
    """Donor-mean hepatocyte scVI latent. Cached to disk after first compute."""
    if HEP_MEAN_CACHE.exists():
        print(f"[hep_scvi] using cache: {HEP_MEAN_CACHE.name}")
        cache = pd.read_csv(HEP_MEAN_CACHE, sep="\t")
        return roster.merge(cache, on=["sample", "dataset"], how="left")

    print(f"[hep_scvi] computing donor-mean latent from {HEP_H5AD.name} (backed)")
    a = ad.read_h5ad(HEP_H5AD, backed="r")
    if "X_scVI" not in a.obsm:
        sys.exit("FATAL: X_scVI missing from hepatocyte atlas; halt.")
    Z = np.asarray(a.obsm["X_scVI"])
    samples = a.obs["sample"].astype(str).values
    datasets = a.obs["dataset"].astype(str).values
    cols = [f"z{i}" for i in range(Z.shape[1])]
    df = pd.DataFrame(Z, columns=cols)
    df["sample"] = samples
    df["dataset"] = datasets
    pb = df.groupby(["sample", "dataset"], observed=True)[cols].mean().reset_index()
    n_hep = (df.groupby(["sample", "dataset"], observed=True)
               .size().rename("n_hepatocytes_for_latent").reset_index())
    pb = pb.merge(n_hep, on=["sample", "dataset"])
    a.file.close()
    pb.to_csv(HEP_MEAN_CACHE, sep="\t", index=False)
    print(f"[hep_scvi] cached {len(pb)} donor rows -> {HEP_MEAN_CACHE.name}")
    return roster.merge(pb, on=["sample", "dataset"], how="left")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ST_DIR.mkdir(parents=True, exist_ok=True)

    roster = pd.read_csv(DONOR_META, sep="\t")[["sample", "dataset"]]
    print(f"[roster] {len(roster)} donors")

    fs = pd.read_csv(FSTAGE_DOC, sep="\t")
    fs["F_stage_documented"] = pd.to_numeric(fs["F_stage_documented"],
                                              errors="coerce")
    labeled = fs.dropna(subset=["F_stage_documented"])[
        ["sample", "dataset", "F_stage_documented"]
    ]
    print(f"[labels] {len(labeled)} donors with documented F-stage")
    print("[labels] class counts:",
          dict(zip(*np.unique(labeled["F_stage_documented"].astype(int),
                              return_counts=True))))

    # ------------------------------------------------------------------
    # Build feature tables on the full roster, then subset to labeled rows
    # for training.
    # ------------------------------------------------------------------
    F_ct = build_celltype_features(roster)
    F_nmf = build_cnmf_features(roster)
    F_pt = build_pseudotime_features(roster)
    F_hep = build_hep_scvi_mean(roster)

    method_results = []
    method_predictions = roster.copy()
    cm_blocks = []

    def fit_eval_logit(name, feat_df, feat_cols, pred_col):
        """Train on labeled subset, LOOCV, predict for everyone with features."""
        merged = feat_df.merge(labeled, on=["sample", "dataset"], how="left")
        # for training, need labels AND non-null features
        train_mask = (merged["F_stage_documented"].notna() &
                      merged[feat_cols].notna().all(axis=1))
        n_train = int(train_mask.sum())
        if n_train < 20:
            print(f"[{name}] only {n_train} usable training donors; skipping")
            method_results.append({
                "method": name,
                "n_train": n_train,
                "qwk_loocv": np.nan,
                "accuracy_loocv": np.nan,
                "off_by_one_loocv": np.nan,
            })
            method_predictions[pred_col] = np.nan
            return

        X_train = merged.loc[train_mask, feat_cols].to_numpy()
        y_train = merged.loc[train_mask, "F_stage_documented"].astype(int).to_numpy()
        print(f"[{name}] training n={n_train}; X shape={X_train.shape}")

        # LOOCV
        y_true, y_pred = loocv_ordinal(X_train, y_train)
        qwk, acc, off1, cm = qwk_metrics(y_true, y_pred)
        print(f"[{name}] LOOCV  QWK={qwk:.3f}  acc={acc:.3f}  off-by-1={off1:.3f}")
        method_results.append({
            "method": name,
            "n_train": n_train,
            "qwk_loocv": round(qwk, 4),
            "accuracy_loocv": round(acc, 4),
            "off_by_one_loocv": round(off1, 4),
        })
        cm_blocks.append((name, cm))

        # refit on all labeled donors, predict where features exist
        clf = make_ordinal_clf()
        clf.fit(X_train, y_train)
        feat_full_mask = merged[feat_cols].notna().all(axis=1)
        X_all = merged.loc[feat_full_mask, feat_cols].to_numpy()
        argmax, _ = predict_with_full_proba(clf, X_all)
        pred_series = pd.Series(np.nan, index=merged.index)
        pred_series.loc[feat_full_mask] = argmax
        # cast to nullable int by passing through float
        method_predictions[pred_col] = pred_series.astype("Float64").astype("Int64")

    def fit_eval_knn(name, feat_df, feat_cols, pred_col, k=5):
        merged = feat_df.merge(labeled, on=["sample", "dataset"], how="left")
        train_mask = (merged["F_stage_documented"].notna() &
                      merged[feat_cols].notna().all(axis=1))
        n_train = int(train_mask.sum())
        if n_train < 20:
            print(f"[{name}] only {n_train} usable training donors; skipping")
            method_results.append({
                "method": name, "n_train": n_train,
                "qwk_loocv": np.nan,
                "accuracy_loocv": np.nan,
                "off_by_one_loocv": np.nan,
            })
            method_predictions[pred_col] = np.nan
            return

        X_train = merged.loc[train_mask, feat_cols].to_numpy()
        y_train = merged.loc[train_mask, "F_stage_documented"].astype(int).to_numpy()
        print(f"[{name}] training n={n_train}; X shape={X_train.shape}")

        y_true, y_pred = loocv_knn_median(X_train, y_train, k=k)
        qwk, acc, off1, cm = qwk_metrics(y_true, y_pred)
        print(f"[{name}] LOOCV  QWK={qwk:.3f}  acc={acc:.3f}  off-by-1={off1:.3f}")
        method_results.append({
            "method": name,
            "n_train": n_train,
            "qwk_loocv": round(qwk, 4),
            "accuracy_loocv": round(acc, 4),
            "off_by_one_loocv": round(off1, 4),
        })
        cm_blocks.append((name, cm))

        feat_full_mask = merged[feat_cols].notna().all(axis=1)
        X_all = merged.loc[feat_full_mask, feat_cols].to_numpy()
        argmax, _ = knn_predict_all(X_train, y_train, X_all, k=k)
        pred_series = pd.Series(np.nan, index=merged.index)
        pred_series.loc[feat_full_mask] = argmax
        method_predictions[pred_col] = pred_series.astype("Float64").astype("Int64")

    # 1. Cell-type proportions
    ct_cols = [f"frac_{lin}" for lin in LINEAGE_MAP]
    fit_eval_logit("celltype_prop", F_ct, ct_cols, "F_stage_celltype")

    # 2. cNMF programs
    nmf_cols = [c for c in F_nmf.columns if c.startswith("cnmf_global_k16_")]
    fit_eval_logit("cnmf_k16", F_nmf, nmf_cols, "F_stage_cnmf")

    # 3. kNN on hep scVI mean
    hep_cols = [c for c in F_hep.columns if c.startswith("z")]
    fit_eval_knn("knn5_scvi_hep", F_hep, hep_cols, "F_stage_knn_scvi", k=5)

    # 4. Pseudotime-only
    pt_cols = ["macrophage_pseudotime_mean",
               "hepatocyte_pseudotime_mean",
               "progressor_frac"]
    fit_eval_logit("pseudotime", F_pt, pt_cols, "F_stage_pseudotime")

    # ------------------------------------------------------------------
    # scVI baseline for the table
    # ------------------------------------------------------------------
    scvi_qwk = np.nan
    if SCVI_PRED.exists():
        s = pd.read_csv(SCVI_PRED, sep="\t")
        if "classifier_qwk_loocv" in s.columns:
            vals = s["classifier_qwk_loocv"].dropna().unique()
            if len(vals) > 0:
                scvi_qwk = float(vals[0])
                # Also bring in argmax predictions
                s_pred = s[["sample", "dataset", "F_stage_predicted_argmax"]]\
                    .rename(columns={"F_stage_predicted_argmax": "F_stage_scvi"})
                method_predictions = method_predictions.merge(
                    s_pred, on=["sample", "dataset"], how="left"
                )
                method_predictions["F_stage_scvi"] = (
                    method_predictions["F_stage_scvi"]
                    .astype("Float64").astype("Int64")
                )
    method_results.append({
        "method": "scvi_ordinal (ref)",
        "n_train": 58,
        "qwk_loocv": scvi_qwk,
        "accuracy_loocv": np.nan,
        "off_by_one_loocv": np.nan,
    })

    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------
    results_df = pd.DataFrame(method_results)
    results_df.to_csv(OUT_COMP, sep="\t", index=False)
    print(f"\n[output] comparison table -> {OUT_COMP}")
    print(results_df.to_string(index=False))

    method_predictions.to_csv(OUT_PRED, sep="\t", index=False)
    print(f"\n[output] alt predictions -> {OUT_PRED}")

    with open(OUT_CMTXT, "w") as f:
        for name, cm in cm_blocks:
            f.write(f"=== {name} ===\n")
            f.write("rows=true F0..F4, cols=pred F0..F4\n")
            f.write(pd.DataFrame(cm,
                                  index=[f"F{i}" for i in range(5)],
                                  columns=[f"F{i}" for i in range(5)]).to_string())
            f.write("\n\n")
    print(f"[output] confusion matrices -> {OUT_CMTXT}")

    # ------------------------------------------------------------------
    # Disagreement summary on the 211 non-Andrews donors (test set)
    # ------------------------------------------------------------------
    test = method_predictions[method_predictions["dataset"] != "GSE202379"].copy()
    pred_cols = [c for c in [
        "F_stage_celltype", "F_stage_cnmf",
        "F_stage_knn_scvi", "F_stage_pseudotime",
        "F_stage_scvi",
    ] if c in test.columns]
    if pred_cols:
        # how often do all 5 methods agree?
        agree_all = test[pred_cols].apply(
            lambda r: r.dropna().nunique() <= 1 and r.notna().sum() == len(pred_cols),
            axis=1,
        )
        # majority vote / spread
        spread = test[pred_cols].apply(
            lambda r: (r.dropna().max() - r.dropna().min())
                       if r.dropna().size > 1 else np.nan,
            axis=1,
        )
        print(f"\n[disagreement] test-set donors n={len(test)} "
              f"(non-GSE202379)")
        print(f"  - donors with full predictions from all {len(pred_cols)} methods: "
              f"{int(test[pred_cols].notna().all(axis=1).sum())}")
        print(f"  - donors where ALL methods agree: {int(agree_all.sum())}")
        print(f"  - max-min spread distribution:\n"
              f"{spread.dropna().astype(int).value_counts().sort_index().to_string()}")


if __name__ == "__main__":
    main()
