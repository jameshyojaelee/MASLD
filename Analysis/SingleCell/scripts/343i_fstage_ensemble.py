#!/usr/bin/env python3
"""343i_fstage_ensemble.py — donor-level ensemble F-stage classifier.

Combines seven pre-computed F-stage predictors at the donor level and reports
LOOCV QWK within the Andrews (GSE202379) labeled set plus external Spearman
rho on non-Andrews donors. Also tracks the Healthy->F4 hallucination rate.

Methods (all donor-level argmax in {0..4}):
  1. scVI ordinal logistic       (donor_fstage_scvi_predicted.tsv)
  2. scVI kNN-5                  (donor_fstage_alt_predictions.tsv : F_stage_knn_scvi)
  3. cNMF k=16                   (donor_fstage_alt_predictions.tsv : F_stage_cnmf)
  4. pseudotime                  (donor_fstage_alt_predictions.tsv : F_stage_pseudotime)
  5. cell-type proportions       (donor_fstage_alt_predictions.tsv : F_stage_celltype)
  6. fibrosis signature          (donor_fstage_fibsignature_predicted.tsv : F_stage_fibsig)
  7. HSC count Rule B            (donor_fstage_hsc_predicted.tsv : F_stage_hsc)

Strategies: mean, median, weighted-mean (LOOCV QWK weights), confidence-weighted
(uses P_max where posterior available; uniform otherwise), stacked ordinal logistic
(LOOCV across Andrews donors).

Outputs:
  donor_fstage_ensemble_predicted.tsv
  ensemble_method_comparison.tsv
  append rows to fstage_method_comparison.tsv
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import cohen_kappa_score

ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
STAGE_DIR = ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
OUT_PRED = STAGE_DIR / "donor_fstage_ensemble_predicted.tsv"
OUT_CMP = STAGE_DIR / "ensemble_method_comparison.tsv"
COMPARE_TSV = STAGE_DIR / "fstage_method_comparison.tsv"

ANDREWS = "GSE202379"

# LOOCV QWK reported in fstage_method_comparison.tsv (used for weighted mean).
METHOD_QWK = {
    "scvi_ordinal":     0.7422,
    "scvi_knn5":        0.7602,
    "cnmf_k16":         0.5628,
    "pseudotime":       0.4135,
    "celltype_prop":    0.3799,
    "fibrosis_signature": 0.3543,
    "hsc_count_ruleB":  0.3014,
}
# Order matters for ensemble stacking
METHODS = list(METHOD_QWK.keys())


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
def load_predictions() -> pd.DataFrame:
    scvi = pd.read_csv(STAGE_DIR / "donor_fstage_scvi_predicted.tsv", sep="\t")
    alt = pd.read_csv(STAGE_DIR / "donor_fstage_alt_predictions.tsv", sep="\t")
    fib = pd.read_csv(STAGE_DIR / "donor_fstage_fibsignature_predicted.tsv", sep="\t")
    hsc = pd.read_csv(STAGE_DIR / "donor_fstage_hsc_predicted.tsv", sep="\t")
    doc = pd.read_csv(STAGE_DIR / "donor_fstage_documented.tsv", sep="\t")
    meta = pd.read_csv(STAGE_DIR / "donor_metadata_extended.tsv", sep="\t")

    keep_meta = ["sample", "dataset", "disease_stage_coarse", "disease_stage_numeric"]
    df = meta[keep_meta].copy()

    df = df.merge(
        scvi[["sample", "F_stage_predicted_argmax",
              "P_F0", "P_F1", "P_F2", "P_F3", "P_F4"]]
            .rename(columns={"F_stage_predicted_argmax": "scvi_ordinal"}),
        on="sample", how="left",
    )

    df = df.merge(
        alt[["sample", "F_stage_knn_scvi", "F_stage_cnmf",
             "F_stage_pseudotime", "F_stage_celltype"]]
            .rename(columns={
                "F_stage_knn_scvi": "scvi_knn5",
                "F_stage_cnmf": "cnmf_k16",
                "F_stage_pseudotime": "pseudotime",
                "F_stage_celltype": "celltype_prop",
            }),
        on="sample", how="left",
    )

    df = df.merge(
        fib[["sample", "F_stage_fibsig",
             "P_F0", "P_F1", "P_F2", "P_F3", "P_F4"]]
            .rename(columns={
                "F_stage_fibsig": "fibrosis_signature",
                "P_F0": "fib_P_F0", "P_F1": "fib_P_F1", "P_F2": "fib_P_F2",
                "P_F3": "fib_P_F3", "P_F4": "fib_P_F4",
            }),
        on="sample", how="left",
    )

    df = df.merge(
        hsc[["sample", "F_stage_hsc"]]
            .rename(columns={"F_stage_hsc": "hsc_count_ruleB"}),
        on="sample", how="left",
    )

    df = df.merge(
        doc[["sample", "F_stage_documented"]],
        on="sample", how="left",
    )

    for m in METHODS:
        df[m] = pd.to_numeric(df[m], errors="coerce")
    df["F_stage_documented"] = pd.to_numeric(df["F_stage_documented"], errors="coerce")
    df["disease_stage_numeric"] = pd.to_numeric(df["disease_stage_numeric"], errors="coerce")
    return df


# ---------------------------------------------------------------------------
# Ensembles
# ---------------------------------------------------------------------------
def _row_scores(row) -> np.ndarray:
    return np.array([row[m] for m in METHODS], dtype=float)


def _clamp(v: float) -> int:
    return int(np.clip(np.rint(v), 0, 4))


def ensemble_mean(df: pd.DataFrame) -> pd.Series:
    return df[METHODS].mean(axis=1, skipna=True).round().clip(0, 4).astype("Int64")


def ensemble_median(df: pd.DataFrame) -> pd.Series:
    return df[METHODS].median(axis=1, skipna=True).round().clip(0, 4).astype("Int64")


def ensemble_weighted_mean(df: pd.DataFrame) -> pd.Series:
    w = np.array([METHOD_QWK[m] if METHOD_QWK[m] >= 0.3 else 0.0 for m in METHODS])
    vals = df[METHODS].to_numpy(dtype=float)
    mask = ~np.isnan(vals)
    W = np.tile(w, (vals.shape[0], 1)) * mask
    num = np.nansum(vals * W, axis=1)
    den = W.sum(axis=1)
    out = np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)
    return pd.Series([_clamp(x) if not np.isnan(x) else pd.NA for x in out],
                     index=df.index, dtype="Int64")


def ensemble_confidence_weighted(df: pd.DataFrame) -> pd.Series:
    """Use per-method P_max when available, else 0.5 (uniform).

    Posterior-bearing methods: scvi_ordinal (P_F0..P_F4) and fibrosis_signature
    (fib_P_F0..fib_P_F4). All others use a constant confidence proxy = their
    LOOCV QWK (i.e., method-level trust), bounded to [0, 1].
    """
    scvi_p = df[["P_F0", "P_F1", "P_F2", "P_F3", "P_F4"]].to_numpy(dtype=float)
    fib_p = df[["fib_P_F0", "fib_P_F1", "fib_P_F2",
                "fib_P_F3", "fib_P_F4"]].to_numpy(dtype=float)
    scvi_pmax = np.nanmax(scvi_p, axis=1)
    fib_pmax = np.nanmax(fib_p, axis=1)

    # Per-method confidence (n_donors x n_methods)
    conf_cols = {}
    for m in METHODS:
        if m == "scvi_ordinal":
            conf_cols[m] = scvi_pmax
        elif m == "fibrosis_signature":
            conf_cols[m] = fib_pmax
        else:
            conf_cols[m] = np.full(len(df), METHOD_QWK[m])

    vals = df[METHODS].to_numpy(dtype=float)
    conf = np.stack([conf_cols[m] for m in METHODS], axis=1)
    mask = ~np.isnan(vals) & ~np.isnan(conf)
    W = np.where(mask, conf, 0.0)
    num = np.nansum(vals * W, axis=1)
    den = W.sum(axis=1)
    out = np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)
    return pd.Series([_clamp(x) if not np.isnan(x) else pd.NA for x in out],
                     index=df.index, dtype="Int64")


def ensemble_stacked(df: pd.DataFrame, andrews_mask: np.ndarray) -> pd.Series:
    """Train a multinomial logistic on Andrews donors via LOOCV; refit on all
    Andrews for predicting non-Andrews donors.

    Missing per-method predictions are imputed with the column mean from the
    training fold (Andrews-only). Output is argmax in {0..4}.
    """
    X_full = df[METHODS].to_numpy(dtype=float)
    y = df["F_stage_documented"].to_numpy(dtype=float)
    n = len(df)

    out = np.full(n, np.nan)

    # ---- LOOCV across Andrews donors (those with documented F-stage) ----
    train_idx = np.where(andrews_mask & ~np.isnan(y))[0]
    y_train_all = y[train_idx].astype(int)

    for i, idx in enumerate(train_idx):
        keep = np.delete(train_idx, i)
        X_tr = X_full[keep].copy()
        y_tr = y[keep].astype(int)
        col_means = np.nanmean(X_tr, axis=0)
        for c in range(X_tr.shape[1]):
            X_tr[np.isnan(X_tr[:, c]), c] = col_means[c]
        X_te = X_full[idx].copy()
        X_te = np.where(np.isnan(X_te), col_means, X_te)
        clf = LogisticRegression(
            solver="lbfgs",
            max_iter=2000, C=1.0,
        )
        # Need >=2 classes; abort fold if degenerate (shouldn't happen).
        if len(np.unique(y_tr)) < 2:
            out[idx] = np.nan
            continue
        clf.fit(X_tr, y_tr)
        pred = clf.predict(X_te.reshape(1, -1))[0]
        out[idx] = int(pred)

    # ---- Refit on full Andrews set and apply to external donors ----
    X_tr_full = X_full[train_idx].copy()
    col_means = np.nanmean(X_tr_full, axis=0)
    for c in range(X_tr_full.shape[1]):
        X_tr_full[np.isnan(X_tr_full[:, c]), c] = col_means[c]
    clf_full = LogisticRegression(
        solver="lbfgs",
        max_iter=2000, C=1.0,
    )
    clf_full.fit(X_tr_full, y_train_all)

    ext_idx = np.where(~andrews_mask)[0]
    X_ext = X_full[ext_idx].copy()
    for c in range(X_ext.shape[1]):
        X_ext[np.isnan(X_ext[:, c]), c] = col_means[c]
    out[ext_idx] = clf_full.predict(X_ext)

    return pd.Series([int(v) if not np.isnan(v) else pd.NA for v in out],
                     index=df.index, dtype="Int64")


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def loocv_qwk_andrews(df: pd.DataFrame, pred_col: str) -> tuple[float, int]:
    """For deterministic ensembles (mean/median/etc.) the 'LOOCV' is just the
    QWK of predictions vs documented labels, since the ensemble itself is not
    fitted. For the stacked ensemble we already populated LOOCV predictions.
    """
    sub = df[(df["dataset"] == ANDREWS) & df["F_stage_documented"].notna()].copy()
    sub = sub[sub[pred_col].notna()]
    if len(sub) < 5:
        return float("nan"), len(sub)
    y_true = sub["F_stage_documented"].astype(int).to_numpy()
    y_pred = sub[pred_col].astype(int).to_numpy()
    return float(cohen_kappa_score(y_true, y_pred, weights="quadratic")), len(sub)


def external_rho(df: pd.DataFrame, pred_col: str) -> tuple[float, float, int]:
    sub = df[(df["dataset"] != ANDREWS) & df["disease_stage_numeric"].notna()].copy()
    sub = sub[sub[pred_col].notna()]
    if len(sub) < 5:
        return float("nan"), float("nan"), len(sub)
    rho, p = spearmanr(sub[pred_col].astype(float),
                       sub["disease_stage_numeric"].astype(float))
    return float(rho), float(p), len(sub)


def healthy_to_f4(df: pd.DataFrame, pred_col: str) -> int:
    """Donors whose disease_stage_coarse == 'Healthy' but predicted F == 4."""
    sub = df[df["disease_stage_coarse"] == "Healthy"].copy()
    sub = sub[sub[pred_col].notna()]
    return int((sub[pred_col].astype(int) == 4).sum())


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    df = load_predictions()
    print(f"[load] {len(df)} donors, "
          f"{(df['dataset'] == ANDREWS).sum()} Andrews, "
          f"{df['F_stage_documented'].notna().sum()} with doc F-stage.")

    # Ensembles
    df["ensemble_mean"] = ensemble_mean(df)
    df["ensemble_median"] = ensemble_median(df)
    df["ensemble_weighted_mean"] = ensemble_weighted_mean(df)
    df["ensemble_confidence_weighted"] = ensemble_confidence_weighted(df)
    andrews_mask = (df["dataset"] == ANDREWS).to_numpy()
    df["ensemble_stacked"] = ensemble_stacked(df, andrews_mask)

    ensembles = [
        "ensemble_mean",
        "ensemble_median",
        "ensemble_weighted_mean",
        "ensemble_confidence_weighted",
        "ensemble_stacked",
    ]

    # Build comparison rows for both ensembles AND individual methods (for context)
    rows = []
    eval_methods = METHODS + ensembles
    for m in eval_methods:
        qwk, n_train = loocv_qwk_andrews(df, m)
        rho, p, n_ext = external_rho(df, m)
        n_h2f4 = healthy_to_f4(df, m)
        rows.append({
            "method": m,
            "n_andrews": n_train,
            "qwk_loocv_andrews": qwk,
            "n_external": n_ext,
            "spearman_rho_external": rho,
            "spearman_p_external": p,
            "healthy_to_f4_count": n_h2f4,
        })
    cmp = pd.DataFrame(rows)
    cmp = cmp.sort_values("qwk_loocv_andrews", ascending=False, na_position="last")

    # --- Save per-donor predictions ---
    keep_cols = (
        ["sample", "dataset", "disease_stage_coarse", "disease_stage_numeric",
         "F_stage_documented"]
        + METHODS
        + ensembles
    )
    df[keep_cols].to_csv(OUT_PRED, sep="\t", index=False)
    print(f"[write] {OUT_PRED}")

    # --- Save ensemble method comparison ---
    cmp.to_csv(OUT_CMP, sep="\t", index=False, float_format="%.4f")
    print(f"[write] {OUT_CMP}")

    # --- Append ensemble rows (only) to fstage_method_comparison.tsv ---
    if COMPARE_TSV.exists():
        existing = pd.read_csv(COMPARE_TSV, sep="\t")
    else:
        existing = pd.DataFrame()
    append_rows = []
    for _, r in cmp.iterrows():
        if r["method"] not in ensembles:
            continue
        # Match the existing schema: method, n_train, qwk_loocv, accuracy_loocv, off_by_one_loocv
        append_rows.append({
            "method": r["method"],
            "n_train": int(r["n_andrews"]) if pd.notna(r["n_andrews"]) else "",
            "qwk_loocv": round(r["qwk_loocv_andrews"], 4)
                if pd.notna(r["qwk_loocv_andrews"]) else "",
            "accuracy_loocv": "",
            "off_by_one_loocv": "",
        })
    if append_rows:
        new_df = pd.DataFrame(append_rows)
        # Replace any existing ensemble rows (idempotent re-run).
        if not existing.empty and "method" in existing.columns:
            existing = existing[~existing["method"].isin(new_df["method"])]
        combined = pd.concat([existing, new_df], ignore_index=True)
        combined.to_csv(COMPARE_TSV, sep="\t", index=False)
        print(f"[update] {COMPARE_TSV} (+{len(new_df)} ensemble rows)")

    # --- Print verdict to stdout ---
    print("\n=== Ensemble vs vanilla scVI ===")
    scvi_qwk = METHOD_QWK["scvi_ordinal"]
    scvi_rho_ref = 0.391  # from task brief
    scvi_h2f4_ref = 27
    for _, r in cmp.iterrows():
        if r["method"] not in ensembles:
            continue
        wins_qwk = (pd.notna(r["qwk_loocv_andrews"])
                    and r["qwk_loocv_andrews"] > scvi_qwk)
        wins_rho = (pd.notna(r["spearman_rho_external"])
                    and r["spearman_rho_external"] > scvi_rho_ref)
        lowers_h = r["healthy_to_f4_count"] < scvi_h2f4_ref
        print(f"  {r['method']:32s}  QWK={r['qwk_loocv_andrews']:.4f}"
              f" rho={r['spearman_rho_external']:.4f}"
              f" H->F4={r['healthy_to_f4_count']:3d}"
              f"   {'QWK+' if wins_qwk else 'QWK-'} "
              f"{'rho+' if wins_rho else 'rho-'} "
              f"{'H4+' if lowers_h else 'H4-'}")


if __name__ == "__main__":
    main()
