#!/usr/bin/env python
"""
343g_hsc_activation_fstage.py

HSC-activation-count-based F-stage classifier for the MASLD scRNA atlas.

Rationale
---------
Histopathology defines fibrosis stage by activated hepatic stellate cell density.
We should be able to recover F-stage from scRNA data by counting how many
fibroblasts/stellate cells per donor are in an "activated myofibroblast" state.
This is fundamentally a counting problem, not a representation-learning problem.

Inputs
------
- Atlas h5ad (backed): scalesc_human_annotated_celltypist.h5ad
- Documented F-stage: donor_fstage_documented.tsv (Andrews 58 donors)
- Donor metadata extended: donor_metadata_extended.tsv

Per-donor features
------------------
- n_fibroblasts, n_activated_fibroblasts, frac_activated_fibroblasts
- frac_activated_normalized = frac_activated / log(n_hepatocytes + 1)
- n_collagen_high_fibroblasts (fraction)
- n_PDGFRB_high_fibroblasts (fraction)

Activation rules (both tried; better LOOCV reported)
----------------------------------------------------
Rule A: mean z(ACTA2,COL1A1,PDGFRB) > 0  AND  mean z(LRAT,GFAP) < 0
Rule B: score_genes(activated_markers) - score_genes(quiescent_markers) > 0

Classifier: ordinal logistic (mord) over per-donor features.
Validation: LOOCV within Andrews 58; external Spearman vs disease_stage_numeric.
Quantile-baseline: frac_activated_fibroblasts alone -> quintile bin -> F0..F4.

Outputs
-------
- donor_fstage_hsc_predicted.tsv
- hsc_method_metrics.tsv
- Appended row in fstage_method_comparison.tsv
"""
from __future__ import annotations
import os
import sys
import gc
from pathlib import Path

import numpy as np
import pandas as pd

import scanpy as sc
import anndata as ad

from scipy import stats
from scipy.sparse import issparse

from sklearn.metrics import cohen_kappa_score, confusion_matrix, accuracy_score
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

try:
    import mord
    HAVE_MORD = True
except ImportError:
    HAVE_MORD = False


PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))

ATLAS_H5AD = PROJECT_ROOT / "Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad"
STAGE_DIR  = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
FSTAGE_DOC = STAGE_DIR / "donor_fstage_documented.tsv"
DONOR_META = STAGE_DIR / "donor_metadata_extended.tsv"
OUT_PRED   = STAGE_DIR / "donor_fstage_hsc_predicted.tsv"
OUT_METRIC = STAGE_DIR / "hsc_method_metrics.tsv"
COMPARE_TSV = STAGE_DIR / "fstage_method_comparison.tsv"

# Marker panels
ACTIVATED_MARKERS = ["ACTA2", "COL1A1", "PDGFRB", "TAGLN", "MYH11"]
QUIESCENT_MARKERS = ["LRAT", "GFAP", "RBP1", "ADAMTSL2", "PPARG"]
COLLAGEN_MARKERS  = ["COL1A1", "COL3A1"]
MIN_FIBROBLASTS_PER_DONOR = 10

# Activation thresholds for collagen/PDGFRB-high (in log1p-normalised space, ~CP10k)
COLLAGEN_THRESHOLD = 1.0   # >=1 in log1p(CP10k) for COL1A1+COL3A1 (per-cell mean)
PDGFRB_THRESHOLD   = 0.5


# ------------------------------------------------------------------
# IO helpers
# ------------------------------------------------------------------

def load_fibroblast_expression():
    """Read fibroblast cells + log-normalised expression for marker genes.

    Returns
    -------
    expr_df : DataFrame indexed by cell barcode; columns = marker gene + sample/dataset
    n_hep_per_donor : Series mapping (sample, dataset) -> n_hepatocytes
    """
    print(f"[load] {ATLAS_H5AD.name} (backed='r')", flush=True)
    a = ad.read_h5ad(str(ATLAS_H5AD), backed="r")
    print(f"[load] shape={a.shape}", flush=True)

    if "cell_type" not in a.obs.columns:
        sys.exit("FATAL: 'cell_type' not in obs columns")
    if "sample" not in a.obs.columns or "dataset" not in a.obs.columns:
        sys.exit("FATAL: 'sample' or 'dataset' not in obs columns")

    ct = a.obs["cell_type"].astype(str)
    fib_mask = ct.eq("Fibroblasts").values
    hep_mask = ct.eq("Hepatocytes").values
    print(f"[load] Fibroblasts={int(fib_mask.sum())}, Hepatocytes={int(hep_mask.sum())}",
          flush=True)

    # n_hepatocytes per donor (used for normalisation feature)
    hep_obs = a.obs.loc[hep_mask, ["sample", "dataset"]].copy()
    n_hep = hep_obs.groupby(["sample", "dataset"], observed=True).size().rename("n_hepatocytes")

    # Materialise fibroblasts to memory
    print("[load] materialising fibroblasts to memory ...", flush=True)
    fib = a[fib_mask].to_memory()
    a.file.close()
    del a
    gc.collect()
    print(f"[load] fibroblast adata: {fib.shape}", flush=True)

    # Use .raw if present (full gene set; pre-HVG)
    if fib.raw is not None:
        print(f"[load] using .raw (n_genes={fib.raw.X.shape[1]})", flush=True)
        score_ad = ad.AnnData(
            X=fib.raw.X.copy(),
            obs=fib.obs.copy(),
            var=fib.raw.var.copy(),
        )
    else:
        print(f"[load] using .X (n_genes={fib.X.shape[1]})", flush=True)
        score_ad = fib

    # Normalise if raw counts
    if issparse(score_ad.X):
        mxv = float(score_ad.X.data.max()) if score_ad.X.nnz > 0 else 0.0
    else:
        mxv = float(score_ad.X.max())
    print(f"[load] max expression value: {mxv:.3f}", flush=True)
    if mxv > 50:
        print("[load] normalize_total + log1p", flush=True)
        sc.pp.normalize_total(score_ad, target_sum=1e4)
        sc.pp.log1p(score_ad)

    # Check marker presence
    gene_pool = set(score_ad.var_names.astype(str))
    all_panels = ACTIVATED_MARKERS + QUIESCENT_MARKERS + COLLAGEN_MARKERS + ["PDGFRB"]
    present = sorted({g for g in all_panels if g in gene_pool})
    missing = sorted({g for g in all_panels if g not in gene_pool})
    print(f"[markers] present: {present}", flush=True)
    print(f"[markers] missing: {missing}", flush=True)

    return score_ad, n_hep


# ------------------------------------------------------------------
# Activation calls
# ------------------------------------------------------------------

def get_gene_matrix(adata, genes):
    """Return dense (n_cells x n_present_genes) matrix for the listed genes."""
    var_names = pd.Index(adata.var_names.astype(str))
    present = [g for g in genes if g in var_names]
    if not present:
        return np.zeros((adata.n_obs, 0)), present
    idx = var_names.get_indexer(present)
    X = adata.X[:, idx]
    if issparse(X):
        X = X.toarray()
    return np.asarray(X, dtype=np.float32), present


def rule_a_calls(adata):
    """Rule A: mean z(activated_core) > 0 AND mean z(quiescent_core) < 0.

    activated_core = {ACTA2, COL1A1, PDGFRB}
    quiescent_core = {LRAT, GFAP}
    z-scores computed across the fibroblast pool.
    """
    act_genes = ["ACTA2", "COL1A1", "PDGFRB"]
    qui_genes = ["LRAT", "GFAP"]
    A, A_present = get_gene_matrix(adata, act_genes)
    Q, Q_present = get_gene_matrix(adata, qui_genes)

    def zscore(M):
        if M.shape[1] == 0:
            return np.zeros((M.shape[0], 0))
        mu = M.mean(axis=0)
        sd = M.std(axis=0)
        sd[sd == 0] = 1.0
        return (M - mu) / sd

    Az = zscore(A); Qz = zscore(Q)
    act_score = Az.mean(axis=1) if Az.shape[1] > 0 else np.zeros(adata.n_obs)
    qui_score = Qz.mean(axis=1) if Qz.shape[1] > 0 else np.zeros(adata.n_obs)
    activated = (act_score > 0) & (qui_score < 0)
    return activated.astype(np.int8), act_score, qui_score, A_present, Q_present


def rule_b_calls(adata):
    """Rule B: sc.tl.score_genes(activated) - sc.tl.score_genes(quiescent) > 0."""
    gene_pool = set(adata.var_names.astype(str))
    act_present = [g for g in ACTIVATED_MARKERS if g in gene_pool]
    qui_present = [g for g in QUIESCENT_MARKERS if g in gene_pool]
    if len(act_present) < 2 or len(qui_present) < 2:
        print("[ruleB] insufficient markers; skipping", flush=True)
        return None, None, None, act_present, qui_present
    sc.tl.score_genes(adata, gene_list=act_present, score_name="_hsc_act_score",
                       use_raw=False, random_state=42)
    sc.tl.score_genes(adata, gene_list=qui_present, score_name="_hsc_qui_score",
                       use_raw=False, random_state=42)
    act_score = adata.obs["_hsc_act_score"].values.astype(np.float32)
    qui_score = adata.obs["_hsc_qui_score"].values.astype(np.float32)
    diff = act_score - qui_score
    activated = (diff > 0).astype(np.int8)
    return activated, act_score, qui_score, act_present, qui_present


def collagen_high_calls(adata):
    M, present = get_gene_matrix(adata, COLLAGEN_MARKERS)
    if M.shape[1] == 0:
        return np.zeros(adata.n_obs, dtype=np.int8), present
    mean_expr = M.mean(axis=1)
    return (mean_expr >= COLLAGEN_THRESHOLD).astype(np.int8), present


def pdgfrb_high_calls(adata):
    M, present = get_gene_matrix(adata, ["PDGFRB"])
    if M.shape[1] == 0:
        return np.zeros(adata.n_obs, dtype=np.int8), present
    return (M[:, 0] >= PDGFRB_THRESHOLD).astype(np.int8), present


# ------------------------------------------------------------------
# Donor-level features
# ------------------------------------------------------------------

def build_donor_features(adata, activated, collagen_high, pdgfrb_high, n_hep):
    """Aggregate per-cell calls into per-donor features."""
    df = pd.DataFrame({
        "sample":  adata.obs["sample"].astype(str).values,
        "dataset": adata.obs["dataset"].astype(str).values,
        "activated":     activated,
        "collagen_high": collagen_high,
        "pdgfrb_high":   pdgfrb_high,
    })
    g = df.groupby(["sample", "dataset"], observed=True)
    feats = pd.DataFrame({
        "n_fibroblasts":             g.size(),
        "n_activated_fibroblasts":   g["activated"].sum(),
        "n_collagen_high_fibroblasts": g["collagen_high"].sum(),
        "n_PDGFRB_high_fibroblasts":   g["pdgfrb_high"].sum(),
    }).reset_index()
    feats["frac_activated_fibroblasts"] = feats["n_activated_fibroblasts"] / feats["n_fibroblasts"]
    feats["frac_collagen_high"] = feats["n_collagen_high_fibroblasts"] / feats["n_fibroblasts"]
    feats["frac_pdgfrb_high"]   = feats["n_PDGFRB_high_fibroblasts"] / feats["n_fibroblasts"]

    # Attach n_hepatocytes
    n_hep_df = n_hep.reset_index()
    feats = feats.merge(n_hep_df, on=["sample", "dataset"], how="left")
    feats["n_hepatocytes"] = feats["n_hepatocytes"].fillna(0).astype(int)
    feats["frac_activated_normalized"] = feats["frac_activated_fibroblasts"] / np.log(feats["n_hepatocytes"] + 1 + 1e-9)
    feats["frac_activated_normalized"] = feats["frac_activated_normalized"].replace([np.inf, -np.inf], np.nan)
    return feats


# ------------------------------------------------------------------
# Classifier helpers
# ------------------------------------------------------------------

FEATURE_COLS = [
    "n_fibroblasts",
    "n_activated_fibroblasts",
    "frac_activated_fibroblasts",
    "frac_activated_normalized",
    "frac_collagen_high",
    "frac_pdgfrb_high",
]


def make_classifier():
    if HAVE_MORD:
        return Pipeline([
            ("scale", StandardScaler()),
            ("clf",   mord.LogisticIT(alpha=1.0)),
        ])
    return Pipeline([
        ("scale", StandardScaler()),
        ("clf",   LogisticRegression(
            solver="lbfgs", C=1.0, max_iter=2000, class_weight="balanced",
        )),
    ])


def loocv_ordinal(X, y):
    n = len(y)
    y_pred = np.empty(n, dtype=int)
    for i in range(n):
        mask = np.ones(n, dtype=bool); mask[i] = False
        clf = make_classifier()
        clf.fit(X[mask], y[mask])
        if HAVE_MORD:
            y_pred[i] = int(clf.predict(X[i:i+1])[0])
        else:
            proba = clf.predict_proba(X[i:i+1])[0]
            y_pred[i] = int(clf.classes_[np.argmax(proba)])
    qwk = cohen_kappa_score(y, y_pred, weights="quadratic", labels=list(range(5)))
    acc = accuracy_score(y, y_pred)
    obo = float(np.mean(np.abs(y - y_pred) <= 1))
    cm = confusion_matrix(y, y_pred, labels=list(range(5)))
    return dict(qwk=float(qwk), accuracy=float(acc), off_by_one=obo,
                y_true=y, y_pred=y_pred, cm=cm)


def quantile_baseline(values, y_true=None):
    """Rank-by-quintile baseline: bin values into 5 quantiles -> F0..F4.

    If y_true given (training/Andrews), report QWK for that subset using the same
    quantile bins computed on the full vector.
    """
    vals = pd.Series(values).astype(float)
    # quantile rank on non-NaN values
    ranks = vals.rank(method="average", pct=True)
    # 5 bins
    bin_id = pd.cut(ranks, bins=[-0.001, 0.2, 0.4, 0.6, 0.8, 1.001],
                     labels=[0, 1, 2, 3, 4]).astype("Int64")
    out = bin_id.to_numpy()  # may contain pd.NA -> object dtype
    if y_true is None:
        return out, None
    # restrict to rows where y_true is present
    mask = pd.Series(y_true).notna().values
    pred_int = np.array([int(x) if pd.notna(x) else -1 for x in out])
    yt = np.asarray(y_true, dtype=float)
    keep = (pred_int >= 0) & mask
    if keep.sum() < 5:
        return out, None
    qwk = cohen_kappa_score(yt[keep].astype(int), pred_int[keep],
                             weights="quadratic", labels=list(range(5)))
    return out, float(qwk)


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def main():
    STAGE_DIR.mkdir(parents=True, exist_ok=True)

    # Load
    adata, n_hep = load_fibroblast_expression()

    # Both rules
    print("\n[ruleA] computing activated calls (z-score rule)", flush=True)
    act_A, A_act_score, A_qui_score, A_actg, A_quig = rule_a_calls(adata)
    print(f"[ruleA] activated cells: {int(act_A.sum())}/{adata.n_obs} "
          f"({100*act_A.mean():.1f}%)", flush=True)

    print("\n[ruleB] computing activated calls (score_genes diff rule)", flush=True)
    res = rule_b_calls(adata)
    act_B, B_act_score, B_qui_score, B_actg, B_quig = res
    if act_B is not None:
        print(f"[ruleB] activated cells: {int(act_B.sum())}/{adata.n_obs} "
              f"({100*act_B.mean():.1f}%)", flush=True)

    # Shared collagen-high / pdgfrb-high calls
    col_high, col_present = collagen_high_calls(adata)
    pdg_high, pdg_present = pdgfrb_high_calls(adata)
    print(f"[derived] collagen-high cells: {int(col_high.sum())} (genes={col_present})",
          flush=True)
    print(f"[derived] PDGFRB-high cells:   {int(pdg_high.sum())} (genes={pdg_present})",
          flush=True)

    # Build donor features for both rules
    feats_A = build_donor_features(adata, act_A, col_high, pdg_high, n_hep)
    feats_A["rule"] = "A"
    if act_B is not None:
        feats_B = build_donor_features(adata, act_B, col_high, pdg_high, n_hep)
        feats_B["rule"] = "B"
    else:
        feats_B = None

    # Attach documented F-stage
    fs = pd.read_csv(FSTAGE_DOC, sep="\t")
    fs["F_stage_documented"] = pd.to_numeric(fs["F_stage_documented"], errors="coerce")
    fs_lab = fs.dropna(subset=["F_stage_documented"])[["sample", "dataset", "F_stage_documented"]]
    fs_lab["F_stage_documented"] = fs_lab["F_stage_documented"].astype(int)
    print(f"\n[label] {len(fs_lab)} documented Andrews donors", flush=True)

    # Gate by min fibroblasts (>=10)
    def train_eval(feats, label):
        f = feats.merge(fs_lab, on=["sample", "dataset"], how="inner")
        f = f[f["n_fibroblasts"] >= MIN_FIBROBLASTS_PER_DONOR].copy()
        f = f.dropna(subset=FEATURE_COLS + ["F_stage_documented"])
        n_train = len(f)
        print(f"  [{label}] n_train (>= {MIN_FIBROBLASTS_PER_DONOR} fibroblasts, "
              f"complete feats): {n_train}", flush=True)
        if n_train < 20:
            return None
        X = f[FEATURE_COLS].to_numpy()
        y = f["F_stage_documented"].astype(int).to_numpy()
        res = loocv_ordinal(X, y)
        res["n_train"] = n_train
        res["train_samples"] = f["sample"].tolist()
        res["train_y"] = y
        return res

    print("\n[loocv] Rule A", flush=True)
    res_A = train_eval(feats_A, "ruleA")
    print(f"  Rule A QWK={res_A['qwk']:.3f}  acc={res_A['accuracy']:.3f}  "
          f"obo={res_A['off_by_one']:.3f}", flush=True)

    res_B = None
    if feats_B is not None:
        print("\n[loocv] Rule B", flush=True)
        res_B = train_eval(feats_B, "ruleB")
        if res_B is not None:
            print(f"  Rule B QWK={res_B['qwk']:.3f}  acc={res_B['accuracy']:.3f}  "
                  f"obo={res_B['off_by_one']:.3f}", flush=True)

    # Pick better rule
    if res_B is None or res_A["qwk"] >= res_B["qwk"]:
        chosen = "A"; chosen_res = res_A; chosen_feats = feats_A
    else:
        chosen = "B"; chosen_res = res_B; chosen_feats = feats_B
    print(f"\n[choose] better rule: {chosen}  (QWK={chosen_res['qwk']:.3f})", flush=True)

    print(f"[loocv] confusion matrix (rows=true F0..F4, cols=pred F0..F4):", flush=True)
    print(chosen_res["cm"], flush=True)

    # Refit on full Andrews using chosen rule -> predict all donors
    train_df = chosen_feats.merge(fs_lab, on=["sample", "dataset"], how="inner")
    train_df = train_df[train_df["n_fibroblasts"] >= MIN_FIBROBLASTS_PER_DONOR].copy()
    train_df = train_df.dropna(subset=FEATURE_COLS + ["F_stage_documented"])
    X_train = train_df[FEATURE_COLS].to_numpy()
    y_train = train_df["F_stage_documented"].astype(int).to_numpy()
    clf = make_classifier()
    clf.fit(X_train, y_train)

    # Predict all donors with sufficient fibroblasts
    pred_df = chosen_feats.copy()
    pred_df["F_stage_hsc"] = pd.Series([pd.NA] * len(pred_df), dtype="Int64")
    elig = (pred_df["n_fibroblasts"] >= MIN_FIBROBLASTS_PER_DONOR) & \
           pred_df[FEATURE_COLS].notna().all(axis=1)
    if elig.any():
        X_all = pred_df.loc[elig, FEATURE_COLS].to_numpy()
        if HAVE_MORD:
            yhat = clf.predict(X_all).astype(int)
        else:
            proba = clf.predict_proba(X_all)
            yhat = clf.classes_[np.argmax(proba, axis=1)].astype(int)
        pred_df.loc[elig, "F_stage_hsc"] = yhat

    # Quantile baseline (sanity check) on frac_activated_fibroblasts.
    # Restrict the quantile binning + evaluation to donors that pass the same
    # gate (>=10 fibroblasts) so it's a like-for-like comparison.
    pred_df = pred_df.reset_index(drop=True)
    elig_q = (pred_df["n_fibroblasts"] >= MIN_FIBROBLASTS_PER_DONOR) & \
              pred_df["frac_activated_fibroblasts"].notna()
    qbin_all = pd.Series([pd.NA] * len(pred_df), dtype="Int64")
    if elig_q.any():
        vals = pred_df.loc[elig_q, "frac_activated_fibroblasts"].astype(float)
        ranks = vals.rank(method="average", pct=True)
        bin_id = pd.cut(ranks, bins=[-0.001, 0.2, 0.4, 0.6, 0.8, 1.001],
                         labels=[0, 1, 2, 3, 4]).astype("Int64")
        qbin_all.loc[elig_q] = bin_id.values
    pred_df["F_stage_hsc_quantile_baseline"] = qbin_all

    # Attach Andrews labels by (sample, dataset)
    pred_with_lab = pred_df.merge(
        fs_lab, on=["sample", "dataset"], how="left",
    )
    is_andrews = pred_with_lab["F_stage_documented"].notna()
    has_qpred  = pred_with_lab["F_stage_hsc_quantile_baseline"].notna()
    mask_q = (is_andrews & has_qpred).values
    if mask_q.sum() >= 5:
        yt = pred_with_lab.loc[mask_q, "F_stage_documented"].astype(int).values
        yp = pred_with_lab.loc[mask_q, "F_stage_hsc_quantile_baseline"].astype(int).values
        qwk_qbin = cohen_kappa_score(yt, yp, weights="quadratic",
                                       labels=list(range(5)))
        acc_qbin = accuracy_score(yt, yp)
        obo_qbin = float(np.mean(np.abs(yt - yp) <= 1))
    else:
        qwk_qbin = float("nan"); acc_qbin = float("nan"); obo_qbin = float("nan")
    print(f"\n[quantile baseline] frac_activated quintile-rank -> F0..F4 on "
          f"{int(mask_q.sum())} Andrews donors: QWK={qwk_qbin:.3f}  "
          f"acc={acc_qbin:.3f}  obo={obo_qbin:.3f}", flush=True)
    pred_df["F_stage_hsc_quantile_baseline"] = pd.Series(qbin_all, dtype="Int64")

    # External validation against disease_stage_numeric
    meta = pd.read_csv(DONOR_META, sep="\t")
    cols_keep = ["sample", "dataset", "disease_stage_coarse", "disease_stage_numeric"]
    cols_keep = [c for c in cols_keep if c in meta.columns]
    pred_ext = pred_df.merge(meta[cols_keep], on=["sample", "dataset"], how="left")

    # Overall Spearman rho (chosen rule prediction vs disease_stage_numeric),
    # restricted to donors NOT in Andrews training set (= external) AND with prediction.
    is_train = pred_ext.merge(fs_lab[["sample", "dataset"]],
                                on=["sample", "dataset"],
                                how="left", indicator=True)["_merge"].eq("both").values
    ext_mask = (~is_train) & pred_ext["F_stage_hsc"].notna() & \
                pred_ext["disease_stage_numeric"].notna()
    n_ext = int(ext_mask.sum())
    if n_ext >= 5:
        rho_ext, p_ext = stats.spearmanr(
            pred_ext.loc[ext_mask, "F_stage_hsc"].astype(int).values,
            pred_ext.loc[ext_mask, "disease_stage_numeric"].astype(float).values,
        )
    else:
        rho_ext, p_ext = float("nan"), float("nan")
    print(f"[external] Spearman rho(F_stage_hsc, disease_stage_numeric) on "
          f"n={n_ext} non-Andrews donors with predictions: rho={rho_ext:.3f}  "
          f"p={p_ext:.2e}", flush=True)

    # Per-dataset rho
    perds_rows = []
    for ds, g in pred_ext.loc[ext_mask].groupby("dataset"):
        if len(g) >= 3 and g["disease_stage_numeric"].nunique() >= 2 and \
           g["F_stage_hsc"].nunique() >= 2:
            r, pv = stats.spearmanr(g["F_stage_hsc"].astype(int),
                                     g["disease_stage_numeric"].astype(float))
        else:
            r, pv = float("nan"), float("nan")
        perds_rows.append(dict(dataset=ds, n=int(len(g)), rho=float(r), p_value=float(pv)))
    perds_df = pd.DataFrame(perds_rows).sort_values("n", ascending=False).reset_index(drop=True)
    print("\n[external] per-dataset Spearman rho:", flush=True)
    print(perds_df.to_string(index=False), flush=True)

    # F4-but-Healthy error rate per dataset
    healthy_F4 = (pred_ext["F_stage_hsc"] == 4) & \
                 (pred_ext["disease_stage_coarse"].astype(str) == "Healthy") & ~is_train
    f4_healthy_per_ds = pred_ext.loc[healthy_F4, "dataset"].value_counts().to_dict()
    n_f4_healthy = int(healthy_F4.sum())
    print(f"\n[errors] F4-predicted-but-Healthy donors (external only): n={n_f4_healthy}",
          flush=True)
    for ds, c in f4_healthy_per_ds.items():
        print(f"   {ds}: {c}", flush=True)

    # ---- Write donor_fstage_hsc_predicted.tsv ----
    out_cols = (
        ["sample", "dataset", "F_stage_hsc", "rule"] +
        FEATURE_COLS +
        ["n_hepatocytes", "F_stage_hsc_quantile_baseline"]
    )
    pred_df["rule"] = chosen
    pred_df[out_cols].to_csv(OUT_PRED, sep="\t", index=False)
    print(f"\n[write] {OUT_PRED}", flush=True)

    # ---- Write hsc_method_metrics.tsv ----
    metric_rows = [
        dict(method="hsc_count_ruleA",
             n_train=res_A["n_train"],
             qwk_loocv=round(res_A["qwk"], 4),
             accuracy_loocv=round(res_A["accuracy"], 4),
             off_by_one_loocv=round(res_A["off_by_one"], 4),
             external_rho=(round(rho_ext, 4) if chosen == "A" else None),
             external_n=(n_ext if chosen == "A" else None),
             n_healthy_F4_misassign=(n_f4_healthy if chosen == "A" else None)),
    ]
    if res_B is not None:
        metric_rows.append(dict(
            method="hsc_count_ruleB",
            n_train=res_B["n_train"],
            qwk_loocv=round(res_B["qwk"], 4),
            accuracy_loocv=round(res_B["accuracy"], 4),
            off_by_one_loocv=round(res_B["off_by_one"], 4),
            external_rho=(round(rho_ext, 4) if chosen == "B" else None),
            external_n=(n_ext if chosen == "B" else None),
            n_healthy_F4_misassign=(n_f4_healthy if chosen == "B" else None)))
    metric_rows.append(dict(
        method="hsc_quantile_baseline",
        n_train=int(mask_q.sum()),
        qwk_loocv=round(qwk_qbin, 4) if not np.isnan(qwk_qbin) else None,
        accuracy_loocv=round(acc_qbin, 4) if not np.isnan(acc_qbin) else None,
        off_by_one_loocv=round(obo_qbin, 4) if not np.isnan(obo_qbin) else None,
        external_rho=None, external_n=None, n_healthy_F4_misassign=None,
    ))
    metrics_df = pd.DataFrame(metric_rows)
    metrics_df.to_csv(OUT_METRIC, sep="\t", index=False)
    print(f"[write] {OUT_METRIC}", flush=True)
    print(metrics_df.to_string(index=False), flush=True)

    # ---- Append row to fstage_method_comparison.tsv ----
    if COMPARE_TSV.exists():
        comp = pd.read_csv(COMPARE_TSV, sep="\t")
    else:
        comp = pd.DataFrame(columns=["method", "n_train", "qwk_loocv",
                                       "accuracy_loocv", "off_by_one_loocv"])
    # remove any prior hsc rows (re-run safety)
    comp = comp[~comp["method"].str.startswith("hsc_", na=False)].copy()
    add_rows = []
    chosen_method = f"hsc_count_rule{chosen}"
    add_rows.append(dict(
        method=chosen_method,
        n_train=chosen_res["n_train"],
        qwk_loocv=round(chosen_res["qwk"], 4),
        accuracy_loocv=round(chosen_res["accuracy"], 4),
        off_by_one_loocv=round(chosen_res["off_by_one"], 4),
    ))
    add_rows.append(dict(
        method="hsc_quantile_baseline",
        n_train=int(mask_q.sum()),
        qwk_loocv=round(qwk_qbin, 4) if not np.isnan(qwk_qbin) else "",
        accuracy_loocv=round(acc_qbin, 4) if not np.isnan(acc_qbin) else "",
        off_by_one_loocv=round(obo_qbin, 4) if not np.isnan(obo_qbin) else "",
    ))
    comp = pd.concat([comp, pd.DataFrame(add_rows)], ignore_index=True)
    comp.to_csv(COMPARE_TSV, sep="\t", index=False)
    print(f"[write] {COMPARE_TSV}", flush=True)
    print(comp.to_string(index=False), flush=True)

    print("\n[done] 343g_hsc_activation_fstage", flush=True)


if __name__ == "__main__":
    main()
