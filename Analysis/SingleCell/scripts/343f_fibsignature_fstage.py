#!/usr/bin/env python
"""
343f_fibsignature_fstage.py

Fibrosis-signature-score-based F-stage classifier for the MASLD scRNA atlas.

Rationale
---------
The scVI-latent ordinal classifier (Script 343b) reaches LOOCV QWK 0.742 within
Andrews (GSE202379) but fails to transfer to most external cohorts (Spearman rho
vs disease_stage_coarse drops to ~0.04 in GSE136103; F4 is hallucinated for 27
healthy donors across Liver_Atlas / GSE185477 / GSE189600). Donor-mean scVI
embeddings appear to reaccumulate dataset batch effects.

This script tests a deliberately crude but cohort-portable alternative:
score every cell against canonical hepatic-stellate / collagen / matrix-remodelling
gene sets (textbook activated-HSC markers), aggregate to per-donor signature
means split by lineage, and train an ordinal classifier on the same 58 Andrews
donors. If the scVI failure mode is batch contamination of the latent space,
canonical marker scores should transfer better, even at the cost of internal
LOOCV QWK.

Inputs
------
- Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad
    Full atlas (1,232,318 cells x 37,533 genes; log1p-normalised X), obs has
    'sample', 'dataset', 'cell_type'.
- Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_documented.tsv
    58 Andrews donors with documented F-stage.
- Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv
    Donor roster + disease_stage_coarse + disease_stage_numeric.

Outputs
-------
- donor_fstage_fibsignature_predicted.tsv
- fibsignature_method_metrics.tsv
- fstage_method_comparison.tsv (appended)
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import cohen_kappa_score, confusion_matrix
from scipy import stats

try:
    import mord  # noqa: F401
    HAVE_MORD = True
except ImportError:
    HAVE_MORD = False

# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))

ATLAS_H5AD = PROJECT_ROOT / "Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad"
ST_DIR     = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
FSTAGE_DOC = ST_DIR / "donor_fstage_documented.tsv"
EXT_META   = ST_DIR / "donor_metadata_extended.tsv"
SCVI_PRED  = ST_DIR / "donor_fstage_scvi_predicted.tsv"

OUT_PRED   = ST_DIR / "donor_fstage_fibsignature_predicted.tsv"
OUT_METRIC = ST_DIR / "fibsignature_method_metrics.tsv"
OUT_COMP   = ST_DIR / "fstage_method_comparison.tsv"
OUT_CONFUS = ST_DIR / "fibsignature_method_confusion.txt"

# Signature gene sets (canonical hepatic fibrosis markers)
SIGNATURES = {
    "fibrosis_core": [
        "ACTA2", "COL1A1", "COL1A2", "COL3A1", "TIMP1", "TIMP2", "LOX", "LOXL2",
        "PDGFRB", "PDGFRA", "FAP", "TGFB1", "SPARC", "VIM", "S100A4", "TAGLN",
        "POSTN", "CTGF",
    ],
    "fibrosis_advanced": [
        "CCN2", "CTHRC1", "INHBA", "MMP2", "SERPINH1", "BGN", "DCN", "LUM",
    ],
    "hsc_activated": [
        "ACTA2", "COL1A1", "COL3A1", "PDGFRB", "RGS5", "MYH11",
    ],
    "hsc_quiescent": [
        "LRAT", "GFAP", "DES", "NGFR", "RBP1", "PPARG", "ADAMTSL2",
    ],
}

# Lineage -> exact strings appearing in obs['cell_type']
LINEAGE_MAP = {
    "Hepatocytes": ["Hepatocytes"],
    "Fibroblasts": ["Fibroblasts"],
    "Macrophages": ["Macrophages", "Mono+mono derived cells"],
    "LSEC":        ["Endothelial cells"],
}

# ---------------------------------------------------------------------------

def make_ordinal_clf():
    if HAVE_MORD:
        import mord
        return Pipeline([("scale", StandardScaler()),
                         ("clf", mord.LogisticIT(alpha=1.0))])
    return Pipeline([("scale", StandardScaler()),
                     ("clf", LogisticRegression(solver="lbfgs",
                                                C=1.0,
                                                max_iter=2000,
                                                class_weight="balanced"))])


def predict_full_proba(clf, X):
    proba = clf.predict_proba(X)
    classes = clf.named_steps["clf"].classes_
    full = np.zeros((proba.shape[0], 5), dtype=float)
    for j, c in enumerate(classes):
        full[:, int(c)] = proba[:, j]
    rs = full.sum(axis=1, keepdims=True)
    rs[rs == 0] = 1.0
    full = full / rs
    return full.argmax(axis=1), full


def loocv(X, y):
    n = len(y)
    y_pred = np.empty(n, dtype=int)
    for i in range(n):
        mask = np.ones(n, dtype=bool); mask[i] = False
        clf = make_ordinal_clf()
        clf.fit(X[mask], y[mask])
        ar, _ = predict_full_proba(clf, X[i:i+1])
        y_pred[i] = int(ar[0])
    return y, y_pred


def score_lineage(adata_full, lineage_labels, ct_col="cell_type"):
    """Subset full atlas (backed) to a lineage in memory and score every signature.

    Returns DataFrame indexed by cell barcode with columns:
       sample, dataset, <one column per signature>.
    """
    mask = adata_full.obs[ct_col].astype(str).isin(lineage_labels).values
    n_cells = int(mask.sum())
    if n_cells == 0:
        return pd.DataFrame(columns=["sample", "dataset"] + list(SIGNATURES.keys()))
    print(f"  [{','.join(lineage_labels)}] n_cells={n_cells:,}; loading into memory...")
    sub = adata_full[mask, :].to_memory()
    # X is already log1p-normalised per task spec.
    sample = sub.obs["sample"].astype(str).values
    dataset = sub.obs["dataset"].astype(str).values

    out = pd.DataFrame({"sample": sample, "dataset": dataset}, index=sub.obs_names)

    for sig_name, glist in SIGNATURES.items():
        present = [g for g in glist if g in sub.var_names]
        if len(present) < 2:
            print(f"    [skip] signature {sig_name}: only {len(present)} markers present")
            out[sig_name] = np.nan
            continue
        # Use scanpy.tl.score_genes (background-corrected, ctrl_size=50)
        sc.tl.score_genes(sub, gene_list=present, score_name=sig_name,
                          ctrl_size=50, n_bins=25, random_state=42,
                          use_raw=False)
        out[sig_name] = sub.obs[sig_name].astype(float).values

    return out


def main():
    OUT_PRED.parent.mkdir(parents=True, exist_ok=True)

    print(f"[load] atlas (backed): {ATLAS_H5AD}")
    a = ad.read_h5ad(ATLAS_H5AD, backed="r")
    print(f"[load] {a.shape[0]:,} cells x {a.shape[1]:,} genes")
    if "cell_type" not in a.obs.columns:
        sys.exit("FATAL: 'cell_type' missing from obs")
    print(f"[load] cell_type distribution:")
    print(a.obs["cell_type"].astype(str).value_counts().head(20).to_string())

    # --- Per-lineage scoring ---
    per_cell_scores = {}
    for lineage, labels in LINEAGE_MAP.items():
        print(f"\n[score] lineage = {lineage}")
        df_cells = score_lineage(a, labels)
        per_cell_scores[lineage] = df_cells

    a.file.close()

    # --- Aggregate to donor x signature x lineage ---
    print("\n[aggregate] mean-score per (sample, dataset, lineage, signature)")
    donor_feat = None
    for lineage, df_cells in per_cell_scores.items():
        if df_cells.empty:
            continue
        agg = (df_cells
               .groupby(["sample", "dataset"], observed=True)[list(SIGNATURES.keys())]
               .mean()
               .add_prefix(f"{lineage}__")
               .reset_index())
        n_cells = (df_cells.groupby(["sample", "dataset"], observed=True)
                   .size().rename(f"n_{lineage}").reset_index())
        agg = agg.merge(n_cells, on=["sample", "dataset"])
        donor_feat = agg if donor_feat is None else donor_feat.merge(agg, on=["sample", "dataset"], how="outer")

    # Derived ratios per task spec.
    # quiescent_to_activated ratio in fibroblasts (mean(hsc_activated) - mean(hsc_quiescent))
    if "Fibroblasts__hsc_activated" in donor_feat.columns:
        donor_feat["fib_activated_minus_quiescent"] = (
            donor_feat["Fibroblasts__hsc_activated"]
            - donor_feat["Fibroblasts__hsc_quiescent"]
        )
    # Stage-progressive ratio across all cells (use weighted mean across lineages).
    # Compute per-cell pooled mean from union of per_cell_scores then group.
    all_cells = pd.concat(list(per_cell_scores.values()), axis=0, ignore_index=False)
    pooled = (all_cells
              .groupby(["sample", "dataset"], observed=True)[["fibrosis_advanced", "hsc_quiescent"]]
              .mean().reset_index())
    pooled["pooled_advanced_minus_quiescent"] = pooled["fibrosis_advanced"] - pooled["hsc_quiescent"]
    donor_feat = donor_feat.merge(
        pooled[["sample", "dataset", "pooled_advanced_minus_quiescent"]],
        on=["sample", "dataset"], how="left",
    )

    # --- Attach labels ---
    fs = pd.read_csv(FSTAGE_DOC, sep="\t")
    fs["F_stage_documented"] = pd.to_numeric(fs["F_stage_documented"], errors="coerce")
    donor_feat = donor_feat.merge(
        fs[["sample", "dataset", "F_stage_documented"]],
        on=["sample", "dataset"], how="left",
    )

    # Feature columns: every signature x lineage mean + 2 ratios.
    feature_cols = [c for c in donor_feat.columns
                    if c not in ("sample", "dataset", "F_stage_documented")
                    and not c.startswith("n_")]
    print(f"[features] n_donor_rows={len(donor_feat)}, n_features={len(feature_cols)}")
    print(f"[features] columns: {feature_cols}")

    # Impute missing features (donors without one of the four lineages get NaN
    # for that lineage's scores). Fill with column mean across labelled donors.
    feat_mat = donor_feat[feature_cols].astype(float).copy()
    label_mask = donor_feat["F_stage_documented"].notna()
    means_train = feat_mat.loc[label_mask].mean(axis=0)
    feat_mat = feat_mat.fillna(means_train)
    donor_feat[feature_cols] = feat_mat

    # --- Training set ---
    train_idx = donor_feat.index[label_mask]
    X_train = donor_feat.loc[train_idx, feature_cols].to_numpy()
    y_train = donor_feat.loc[train_idx, "F_stage_documented"].astype(int).to_numpy()
    n_train = len(y_train)
    print(f"\n[train] n_train={n_train}  X shape={X_train.shape}")
    print(f"[train] class counts: {dict(zip(*np.unique(y_train, return_counts=True)))}")

    # --- LOOCV ---
    print("[loocv] leave-one-donor-out...")
    y_true, y_pred = loocv(X_train, y_train)
    qwk = float(cohen_kappa_score(y_true, y_pred, weights="quadratic",
                                  labels=list(range(5))))
    acc = float(np.mean(y_true == y_pred))
    off1 = float(np.mean(np.abs(y_true - y_pred) <= 1))
    cm = confusion_matrix(y_true, y_pred, labels=list(range(5)))
    print(f"[loocv] QWK={qwk:.4f}  acc={acc:.4f}  off-by-one={off1:.4f}")
    print("[loocv] confusion (rows=true F0..F4, cols=pred F0..F4):")
    print(cm)

    # --- Refit on all 58 train donors; predict for everyone ---
    clf = make_ordinal_clf()
    clf.fit(X_train, y_train)
    X_all = donor_feat[feature_cols].to_numpy()
    pred_argmax, pred_proba = predict_full_proba(clf, X_all)

    out = donor_feat[["sample", "dataset"] + feature_cols + ["F_stage_documented"]].copy()
    out["F_stage_fibsig"] = pred_argmax.astype(int)
    for k in range(5):
        out[f"P_F{k}"] = pred_proba[:, k].round(6)
    out["classifier_qwk_loocv"] = round(qwk, 4)

    out.to_csv(OUT_PRED, sep="\t", index=False)
    print(f"[write] {OUT_PRED}  ({len(out)} rows)")

    # --- External validation against extended meta ---
    print("\n[external] merging with donor_metadata_extended for external rho")
    meta = pd.read_csv(EXT_META, sep="\t")
    merged = out.merge(
        meta[["sample", "dataset", "disease_stage_coarse", "disease_stage_numeric", "F_stage_doc_source"]],
        on=["sample", "dataset"], how="left",
    )

    # external = donors WITHOUT documented F-stage (i.e. not Andrews training)
    external_mask = merged["F_stage_documented"].isna() & merged["disease_stage_numeric"].notna()
    ext = merged.loc[external_mask].copy()
    n_ext = len(ext)
    if n_ext >= 3:
        rho_ext, p_ext = stats.spearmanr(
            ext["F_stage_fibsig"].astype(float),
            pd.to_numeric(ext["disease_stage_numeric"], errors="coerce"),
        )
        rho_ext = float(rho_ext); p_ext = float(p_ext)
    else:
        rho_ext, p_ext = np.nan, np.nan

    # Healthy_F4 misassignment count
    healthy_mask = (merged["disease_stage_coarse"] == "Healthy") & merged["F_stage_documented"].isna()
    n_healthy_f4 = int(((merged.loc[healthy_mask, "F_stage_fibsig"]) == 4).sum())
    n_healthy_total = int(healthy_mask.sum())

    # Per-dataset rho
    print(f"[external] n_external_donors_with_stage = {n_ext}")
    print(f"[external] overall Spearman rho(F_stage_fibsig, disease_stage_numeric) = {rho_ext:.4f}  (p={p_ext:.2e})")
    print(f"[external] healthy donors predicted F4: {n_healthy_f4} / {n_healthy_total}")
    perds_rows = []
    for ds, g in ext.groupby("dataset"):
        if g["F_stage_fibsig"].nunique() < 2 or g["disease_stage_numeric"].nunique() < 2:
            r, p = np.nan, np.nan
        else:
            r, p = stats.spearmanr(g["F_stage_fibsig"].astype(float),
                                    pd.to_numeric(g["disease_stage_numeric"], errors="coerce"))
        perds_rows.append(dict(dataset=ds, n_donors=int(len(g)),
                               rho=float(r) if not np.isnan(r) else np.nan,
                               p_value=float(p) if not np.isnan(p) else np.nan,
                               n_F4_predicted=int((g["F_stage_fibsig"]==4).sum()),
                               n_Healthy=int((g["disease_stage_coarse"]=="Healthy").sum())))
    perds = pd.DataFrame(perds_rows).sort_values("n_donors", ascending=False)

    # Cross-tab predicted vs disease_stage_coarse (external donors)
    if n_ext > 0:
        crosstab_ext = pd.crosstab(ext["F_stage_fibsig"], ext["disease_stage_coarse"]).reindex(
            index=list(range(5)),
            columns=["Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"],
            fill_value=0,
        )
    else:
        crosstab_ext = pd.DataFrame()

    # --- Comparison row ---
    metric_row = pd.DataFrame([dict(
        method="fibrosis_signature",
        n_train=n_train,
        qwk_loocv=round(qwk, 4),
        accuracy=round(acc, 4),
        off_by_one=round(off1, 4),
        external_rho_disease_stage_numeric=round(rho_ext, 4) if not np.isnan(rho_ext) else np.nan,
        external_pvalue=round(p_ext, 6) if not np.isnan(p_ext) else np.nan,
        n_external_donors=int(n_ext),
        n_healthy_F4_misassign=n_healthy_f4,
        n_healthy_total=n_healthy_total,
    )])
    metric_row.to_csv(OUT_METRIC, sep="\t", index=False)
    print(f"[write] {OUT_METRIC}")

    # Append to method_comparison.tsv (preserving existing rows; harmonise schema)
    if OUT_COMP.exists():
        existing = pd.read_csv(OUT_COMP, sep="\t")
    else:
        existing = pd.DataFrame()
    add_row = pd.DataFrame([dict(
        method="fibrosis_signature",
        n_train=n_train,
        qwk_loocv=round(qwk, 4),
        accuracy_loocv=round(acc, 4),
        off_by_one_loocv=round(off1, 4),
    )])
    # drop any prior fibrosis_signature row to keep idempotent
    if "method" in existing.columns:
        existing = existing[existing["method"] != "fibrosis_signature"]
    new_comp = pd.concat([existing, add_row], axis=0, ignore_index=True, sort=False)
    new_comp.to_csv(OUT_COMP, sep="\t", index=False)
    print(f"[write] {OUT_COMP} (now {len(new_comp)} rows)")

    # --- Write confusion / external summary as text ---
    lines = []
    L = lines.append
    L("# fibsignature F-stage classifier confusion / external summary")
    L(f"# n_train={n_train}; n_features={len(feature_cols)}")
    L(f"LOOCV QWK = {qwk:.4f}  acc = {acc:.4f}  off-by-one = {off1:.4f}")
    L("LOOCV confusion (rows=true F0..F4, cols=pred F0..F4):")
    L(str(cm))
    L("")
    L(f"External n={n_ext}; overall Spearman rho={rho_ext:.4f} (p={p_ext:.2e})")
    L(f"Healthy donors predicted F4: {n_healthy_f4} / {n_healthy_total}")
    L("Per-dataset (external):")
    L(perds.to_string(index=False))
    L("")
    L("External cross-tab F_stage_fibsig vs disease_stage_coarse:")
    L(crosstab_ext.to_string())
    OUT_CONFUS.write_text("\n".join(lines) + "\n")
    print(f"[write] {OUT_CONFUS}")

    # --- Compare to existing scVI methods for context ---
    if SCVI_PRED.exists():
        scv = pd.read_csv(SCVI_PRED, sep="\t")
        scv_merge = scv.merge(meta[["sample", "dataset", "disease_stage_numeric", "F_stage_doc_source"]],
                              on=["sample", "dataset"], how="left")
        ext_scv = scv_merge[scv_merge["F_stage_doc_source"].isna() | (scv_merge["F_stage_doc_source"].astype(str)=="")]
        ext_scv = ext_scv.dropna(subset=["F_stage_predicted_argmax", "disease_stage_numeric"])
        if len(ext_scv) >= 3:
            rho_scv, p_scv = stats.spearmanr(
                ext_scv["F_stage_predicted_argmax"].astype(float),
                pd.to_numeric(ext_scv["disease_stage_numeric"], errors="coerce"),
            )
            print(f"[compare] scVI logistic external rho (n={len(ext_scv)}): {rho_scv:.4f} p={p_scv:.2e}")
            print(f"[compare] fibsig  external rho (n={n_ext}): {rho_ext:.4f} p={p_ext:.2e}")

    print("\n[done] 343f fibsignature F-stage classifier")


if __name__ == "__main__":
    main()
