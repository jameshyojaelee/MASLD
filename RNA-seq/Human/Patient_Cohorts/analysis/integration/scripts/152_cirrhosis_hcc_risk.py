# WARNING: Cirrhosis risk head produces AUROC 0.431 (below chance) on aggregated LOCO predictions. Circular pseudo-labels. HCC molecular score portion is valid.
#!/usr/bin/env python3
"""
152_cirrhosis_hcc_risk.py — Cirrhosis + HCC Molecular Risk
-------------------------------------------------------------
Two complementary risk heads built on the prognosis feature matrix:

1. **Cirrhosis Risk Head**: Elastic net logistic regression predicting
   binarized P(F4) among F2-F3 samples (the clinical "at risk" group).
   Validated with LOCO-CV (6 fibrosis folds).

2. **HCC Molecular Predisposition Head**: Gene-set score derived from
   HCC COLOC genes (PP.H4-weighted z-scored expression). NOT a classifier—
   this quantifies germline HCC risk reflected in liver transcriptome.
   Validated against fibrosis stage, P(F4), and NMF subtype.

LIMITATIONS (honest framing):
  - Cirrhosis risk model uses pseudo-labels (P(F4) from optimal transport),
    NOT actual clinical progression outcomes.
  - HCC score reflects transcriptomic expression of genetically causal HCC
    loci, not prospective HCC risk. Cross-sectional design cannot establish
    temporal precedence.
  - Both scores require external longitudinal validation.

Input files (from Script 150):
  - results/prognosis/prognosis_feature_matrix.csv
  - results/prognosis/prognosis_pseudo_labels.csv
  - results/prognosis/hcc_coloc_genes.csv
  - results/staging_classifier/modeling_metadata.csv
  - RNA-seq/results/subtypes/nmf_assignments.csv

Output (to results/prognosis/):
  - cirrhosis_risk_loco.csv       (per-sample LOCO predictions for F2-F3)
  - hcc_molecular_score.csv       (per-sample HCC molecular score)
  - hcc_validation.csv            (HCC score vs fibrosis/P(F4)/subtype)

SLURM: io, 4 CPUs, 8G RAM, 48h
Env:   micromamba activate spatial
"""

import logging
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INTEG = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RNASEQ = os.path.join(BASE, "RNA-seq")

FEATURE_MATRIX = os.path.join(INTEG, "results/prognosis/prognosis_feature_matrix.csv")
PSEUDO_LABELS = os.path.join(INTEG, "results/prognosis/prognosis_pseudo_labels.csv")
HCC_COLOC_GENES = os.path.join(INTEG, "results/prognosis/hcc_coloc_genes.csv")
MODELING_META = os.path.join(INTEG, "results/staging_classifier/modeling_metadata.csv")
NMF_ASSIGNMENTS = os.path.join(RNASEQ, "results/subtypes/nmf_assignments.csv")

OUT_DIR = os.path.join(INTEG, "results/prognosis")


# ── Utility ──────────────────────────────────────────────────────────────────
def safe_load(path, **kwargs):
    """Load CSV with existence check."""
    if not os.path.exists(path):
        log.warning("File not found: %s", path)
        return None
    log.info("Loading %s", os.path.basename(path))
    return pd.read_csv(path, **kwargs)


# ══════════════════════════════════════════════════════════════════════════════
# PART 1: Cirrhosis Risk Head (Elastic Net on F2-F3 samples)
# ══════════════════════════════════════════════════════════════════════════════
def cirrhosis_risk_model(features, labels, meta):
    """
    Elastic net logistic regression predicting binarized P(F4) among F2-F3.

    Target population: F2-F3 samples (clinical "gray zone" for cirrhosis risk).
    Pseudo-label: P(F4) binarized at F2-F3 median.
    Validation: LOCO-CV using existing fibrosis folds.
    """
    log.info("=" * 70)
    log.info("PART 1: Cirrhosis Risk Head (F2-F3 samples)")
    log.info("=" * 70)

    # Merge features + labels + metadata
    df = features.join(labels[["p_f4", "fibrosis_stage", "loco_fold_fibrosis"]])

    # Filter to F2-F3 samples with valid P(F4)
    fib_stage = pd.to_numeric(df["fibrosis_stage"], errors="coerce")
    mask_f2f3 = fib_stage.isin([2, 3])
    mask_pf4 = df["p_f4"].notna()
    mask_fold = (df["loco_fold_fibrosis"] != "excluded") & df["loco_fold_fibrosis"].notna()
    mask = mask_f2f3 & mask_pf4 & mask_fold
    df_sub = df.loc[mask].copy()
    n_f2f3 = mask_f2f3.sum()
    n_valid = len(df_sub)
    log.info("F2-F3 samples: %d total, %d with valid P(F4) + LOCO fold", n_f2f3, n_valid)

    if n_valid < 30:
        log.warning("Too few F2-F3 samples (%d < 30). Skipping cirrhosis risk model.", n_valid)
        return pd.DataFrame()

    # Binarize P(F4) at F2-F3 median
    pf4_median = df_sub["p_f4"].median()
    df_sub["cirrhosis_risk_label"] = (df_sub["p_f4"] >= pf4_median).astype(int)
    log.info("P(F4) binarization threshold (median): %.4f", pf4_median)
    log.info("Label distribution: high_risk=%d, low_risk=%d",
             df_sub["cirrhosis_risk_label"].sum(),
             (1 - df_sub["cirrhosis_risk_label"]).sum())

    # Select features: focused set for this model
    # T5 transition score F3->F4, T2 stellate fraction + binary, T1 div genes, T6 HCC score
    priority_features = []
    # Transition scores
    for col in ["trans_F3_to_F4", "trans_F2_to_F3"]:
        if col in features.columns:
            priority_features.append(col)
    # Stellate features
    for col in ["ct_Stellate", "ct_stellate_gt4pct", "ct_stellate_hepatocyte_ratio",
                "ct_fibrogenic_axis", "ct_immune_infiltration"]:
        if col in features.columns:
            priority_features.append(col)
    # HCC genetic score
    if "hcc_genetic_score" in features.columns:
        priority_features.append("hcc_genetic_score")
    # GAS6::MERTK
    for col in ["gas6_mertk_product", "gas6_mertk_ratio"]:
        if col in features.columns:
            priority_features.append(col)
    # Top divergence genes (first 20 by column order = by |d| rank)
    div_cols = [c for c in features.columns if c.startswith("div_")][:20]
    priority_features.extend(div_cols)
    # Clinical sex
    if "clin_sex" in features.columns:
        priority_features.append("clin_sex")

    # Remove any features not present in df_sub
    feature_cols = [c for c in priority_features if c in df_sub.columns]
    log.info("Using %d features for cirrhosis risk model", len(feature_cols))

    X = df_sub[feature_cols].values.astype(float)
    y = df_sub["cirrhosis_risk_label"].values
    folds = df_sub["loco_fold_fibrosis"].values

    # Handle NaN in features: impute with column median
    for j in range(X.shape[1]):
        col_vals = X[:, j]
        nan_mask = np.isnan(col_vals)
        if nan_mask.any():
            median_val = np.nanmedian(col_vals)
            X[nan_mask, j] = median_val if not np.isnan(median_val) else 0.0

    # LOCO-CV
    unique_folds = sorted(set(folds))
    log.info("LOCO folds: %s", unique_folds)

    results = []
    fold_metrics = []

    for fold in unique_folds:
        test_mask = folds == fold
        train_mask = ~test_mask

        if test_mask.sum() == 0 or train_mask.sum() == 0:
            continue

        X_train, X_test = X[train_mask], X[test_mask]
        y_train, y_test = y[train_mask], y[test_mask]

        # Check class balance in train
        if len(np.unique(y_train)) < 2:
            log.warning("Fold %s: only one class in training. Skipping.", fold)
            continue

        # Standardize
        scaler = StandardScaler()
        X_train_s = scaler.fit_transform(X_train)
        X_test_s = scaler.transform(X_test)

        # Elastic net logistic regression (l1_ratio=0.5 for elastic net)
        model = LogisticRegressionCV(
            penalty="elasticnet",
            solver="saga",
            l1_ratios=[0.1, 0.3, 0.5, 0.7, 0.9],
            Cs=10,
            cv=3,
            max_iter=5000,
            random_state=42,
            scoring="roc_auc",
            n_jobs=-1,
        )
        model.fit(X_train_s, y_train)
        proba = model.predict_proba(X_test_s)[:, 1]

        # Per-fold metrics
        fold_auroc = np.nan
        fold_auprc = np.nan
        if len(np.unique(y_test)) == 2:
            fold_auroc = roc_auc_score(y_test, proba)
            fold_auprc = average_precision_score(y_test, proba)
        fold_metrics.append({
            "fold": fold,
            "n_test": test_mask.sum(),
            "auroc": fold_auroc,
            "auprc": fold_auprc,
            "n_nonzero_coefs": np.sum(model.coef_[0] != 0),
        })
        log.info("Fold %s: n=%d, AUROC=%.3f, AUPRC=%.3f, nonzero=%d",
                 fold, test_mask.sum(), fold_auroc, fold_auprc,
                 np.sum(model.coef_[0] != 0))

        # Store per-sample predictions
        test_samples = df_sub.index[test_mask]
        for i, sid in enumerate(test_samples):
            results.append({
                "sample_id": sid,
                "loco_fold": fold,
                "cirrhosis_risk_prob": proba[i],
                "cirrhosis_risk_label": y_test[i],
                "fibrosis_stage": fib_stage.loc[sid],
                "p_f4": df_sub.loc[sid, "p_f4"],
            })

    if len(results) == 0:
        log.warning("No valid LOCO results produced.")
        return pd.DataFrame()

    results_df = pd.DataFrame(results)

    # Aggregate metrics
    metrics_df = pd.DataFrame(fold_metrics)
    mean_auroc = metrics_df["auroc"].mean()
    mean_auprc = metrics_df["auprc"].mean()
    log.info("-" * 50)
    log.info("LOCO Cirrhosis Risk: mean AUROC=%.3f, mean AUPRC=%.3f",
             mean_auroc, mean_auprc)
    log.info("Per-fold breakdown:\n%s", metrics_df.to_string(index=False))

    # Fit final model on all data for feature importances
    scaler_full = StandardScaler()
    X_full_s = scaler_full.fit_transform(X)
    model_full = LogisticRegressionCV(
        penalty="elasticnet",
        solver="saga",
        l1_ratios=[0.1, 0.3, 0.5, 0.7, 0.9],
        Cs=10,
        cv=3,
        max_iter=5000,
        random_state=42,
        scoring="roc_auc",
        n_jobs=-1,
    )
    model_full.fit(X_full_s, y)
    coefs = model_full.coef_[0]
    coef_df = pd.DataFrame({
        "feature": feature_cols,
        "coefficient": coefs,
        "abs_coefficient": np.abs(coefs),
    }).sort_values("abs_coefficient", ascending=False)
    nonzero_coefs = coef_df[coef_df["coefficient"] != 0]
    log.info("Full model: %d / %d nonzero features", len(nonzero_coefs), len(feature_cols))
    log.info("Top 10 features:\n%s", nonzero_coefs.head(10).to_string(index=False))

    # Add metrics row as metadata in output
    results_df.attrs["mean_auroc"] = mean_auroc
    results_df.attrs["mean_auprc"] = mean_auprc

    return results_df


# ══════════════════════════════════════════════════════════════════════════════
# PART 2: HCC Molecular Predisposition Score
# ══════════════════════════════════════════════════════════════════════════════
def hcc_molecular_score(features, labels, meta):
    """
    HCC molecular predisposition score from COLOC genes.

    NOT a classifier. This is a gene-set score:
      score = sum(z-scored expression * PP.H4 weight) for HCC COLOC genes

    The HCC genetic score was already computed in Script 150 (T6).
    Here we load it, validate against clinical annotations, and save.
    """
    log.info("=" * 70)
    log.info("PART 2: HCC Molecular Predisposition Score")
    log.info("=" * 70)

    # Load HCC COLOC genes
    hcc_genes = safe_load(HCC_COLOC_GENES)
    if hcc_genes is None or len(hcc_genes) == 0:
        log.warning("No HCC COLOC genes found. Skipping HCC score.")
        return pd.DataFrame(), pd.DataFrame()

    n_total = len(hcc_genes)
    n_mapped = hcc_genes["in_expression_matrix"].sum()
    log.info("HCC COLOC genes: %d total, %d mapped to expression", n_total, n_mapped)

    # HCC score already in feature matrix (hcc_genetic_score)
    if "hcc_genetic_score" not in features.columns:
        log.warning("hcc_genetic_score not in feature matrix. Skipping.")
        return pd.DataFrame(), pd.DataFrame()

    hcc_score_series = features["hcc_genetic_score"].copy()
    # Reset to numpy for consistent indexing in validation
    hcc_score = hcc_score_series.values.astype(float)
    log.info("HCC score: mean=%.4f, std=%.4f, %d non-NaN",
             np.nanmean(hcc_score), np.nanstd(hcc_score),
             (~np.isnan(hcc_score)).sum())

    sample_ids = features.index.values

    # Build per-sample output
    score_df = pd.DataFrame({
        "sample_id": sample_ids,
        "hcc_molecular_score": hcc_score,
    })
    # Add metadata from labels
    score_df["fibrosis_stage"] = labels["fibrosis_stage"].values
    score_df["nas_score"] = labels["nas_score"].values
    score_df["dataset"] = labels["dataset"].values
    score_df["s2_binary"] = labels["s2_binary"].values
    score_df["p_f4"] = labels["p_f4"].values
    score_df["cps"] = labels["cps"].values

    # Load NMF assignments for subtype annotation
    nmf = safe_load(NMF_ASSIGNMENTS)
    if nmf is not None:
        nmf_map = nmf.set_index("sample_id")["nmf_subtype"]
        score_df["nmf_subtype"] = score_df["sample_id"].map(nmf_map)
    else:
        score_df["nmf_subtype"] = np.nan

    # ── Validation ────────────────────────────────────────────────────────
    log.info("--- HCC Score Validation ---")
    validation_rows = []

    # All validation uses numpy arrays (consistent indexing, no pandas alignment issues)
    fib_arr = pd.to_numeric(score_df["fibrosis_stage"], errors="coerce").values
    pf4_arr = pd.to_numeric(score_df["p_f4"], errors="coerce").values
    cps_arr = pd.to_numeric(score_df["cps"], errors="coerce").values
    nas_arr = pd.to_numeric(score_df["nas_score"], errors="coerce").values
    subtype_arr = score_df["nmf_subtype"].values if "nmf_subtype" in score_df.columns else None

    # 1. Correlation with fibrosis stage
    valid = ~np.isnan(hcc_score) & ~np.isnan(fib_arr) & (fib_arr >= 0)
    if valid.sum() > 10:
        rho, pval = stats.spearmanr(hcc_score[valid], fib_arr[valid])
        validation_rows.append({
            "test": "HCC_score_vs_fibrosis_stage",
            "method": "Spearman",
            "statistic": rho,
            "p_value": pval,
            "n_samples": int(valid.sum()),
            "interpretation": "positive = HCC risk increases with fibrosis",
        })
        log.info("HCC score vs fibrosis stage: rho=%.3f, p=%.2e, n=%d", rho, pval, valid.sum())

    # 2. Correlation with P(F4)
    valid = ~np.isnan(hcc_score) & ~np.isnan(pf4_arr)
    if valid.sum() > 10:
        rho, pval = stats.spearmanr(hcc_score[valid], pf4_arr[valid])
        validation_rows.append({
            "test": "HCC_score_vs_P_F4",
            "method": "Spearman",
            "statistic": rho,
            "p_value": pval,
            "n_samples": int(valid.sum()),
            "interpretation": "positive = HCC risk co-tracks with cirrhosis fate probability",
        })
        log.info("HCC score vs P(F4): rho=%.3f, p=%.2e, n=%d", rho, pval, valid.sum())

    # 3. Correlation with CPS
    valid = ~np.isnan(hcc_score) & ~np.isnan(cps_arr)
    if valid.sum() > 10:
        rho, pval = stats.spearmanr(hcc_score[valid], cps_arr[valid])
        validation_rows.append({
            "test": "HCC_score_vs_CPS",
            "method": "Spearman",
            "statistic": rho,
            "p_value": pval,
            "n_samples": int(valid.sum()),
            "interpretation": "positive = HCC risk tracks with composite progression",
        })
        log.info("HCC score vs CPS: rho=%.3f, p=%.2e, n=%d", rho, pval, valid.sum())

    # 4. Kruskal-Wallis across fibrosis stages
    groups = []
    group_labels = []
    fib_valid_mask = ~np.isnan(fib_arr) & (fib_arr >= 0) & ~np.isnan(hcc_score)
    for stage in sorted(np.unique(fib_arr[~np.isnan(fib_arr) & (fib_arr >= 0)])):
        mask = (fib_arr == stage) & ~np.isnan(hcc_score)
        if mask.sum() > 0:
            groups.append(hcc_score[mask])
            group_labels.append(f"F{int(stage)}")
    if len(groups) >= 2:
        h_stat, h_pval = stats.kruskal(*groups)
        validation_rows.append({
            "test": "HCC_score_Kruskal_Wallis_by_fibrosis",
            "method": "Kruskal-Wallis",
            "statistic": h_stat,
            "p_value": h_pval,
            "n_samples": sum(len(g) for g in groups),
            "interpretation": f"Groups: {group_labels}; significant = HCC score differs by fibrosis",
        })
        log.info("KW across fibrosis stages: H=%.2f, p=%.2e, groups=%s",
                 h_stat, h_pval, group_labels)
        # Per-group means
        for label, group in zip(group_labels, groups):
            log.info("  %s: mean=%.3f, median=%.3f, n=%d",
                     label, np.mean(group), np.median(group), len(group))

    # 5. S2 subtype comparison (Mann-Whitney)
    if subtype_arr is not None:
        s1_mask = (subtype_arr == "S1") & ~np.isnan(hcc_score)
        s2_mask = (subtype_arr == "S2") & ~np.isnan(hcc_score)
        if s1_mask.sum() > 5 and s2_mask.sum() > 5:
            s1_scores = hcc_score[s1_mask]
            s2_scores = hcc_score[s2_mask]
            u_stat, u_pval = stats.mannwhitneyu(s2_scores, s1_scores, alternative="two-sided")
            validation_rows.append({
                "test": "HCC_score_S2_vs_S1",
                "method": "Mann-Whitney U",
                "statistic": u_stat,
                "p_value": u_pval,
                "n_samples": int(s1_mask.sum() + s2_mask.sum()),
                "interpretation": (
                    f"S2 mean={np.mean(s2_scores):.3f} vs S1 mean={np.mean(s1_scores):.3f}; "
                    f"S2 is female-enriched progressor subtype"
                ),
            })
            log.info("HCC score S2 vs S1: U=%.1f, p=%.2e (S2 mean=%.3f, S1 mean=%.3f)",
                     u_stat, u_pval, np.mean(s2_scores), np.mean(s1_scores))

    # 6. NAS score correlation
    valid = ~np.isnan(hcc_score) & ~np.isnan(nas_arr)
    if valid.sum() > 10:
        rho, pval = stats.spearmanr(hcc_score[valid], nas_arr[valid])
        validation_rows.append({
            "test": "HCC_score_vs_NAS",
            "method": "Spearman",
            "statistic": rho,
            "p_value": pval,
            "n_samples": int(valid.sum()),
            "interpretation": "HCC risk vs NAFLD Activity Score",
        })
        log.info("HCC score vs NAS: rho=%.3f, p=%.2e, n=%d", rho, pval, valid.sum())

    validation_df = pd.DataFrame(validation_rows)
    return score_df, validation_df


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════
def main():
    t0 = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    log.info("=" * 70)
    log.info("Script 152: Cirrhosis + HCC Molecular Risk")
    log.info("=" * 70)
    log.info("")
    log.info("IMPORTANT CAVEATS:")
    log.info("  - Cirrhosis risk uses pseudo-labels from optimal transport, NOT clinical outcomes.")
    log.info("  - HCC score reflects transcriptomic expression of HCC COLOC loci,")
    log.info("    NOT prospective HCC risk. Cross-sectional design only.")
    log.info("  - Both scores require external longitudinal validation.")
    log.info("")

    # Load data
    features = safe_load(FEATURE_MATRIX, index_col="sample_id")
    labels = safe_load(PSEUDO_LABELS, index_col="sample_id")
    meta = safe_load(MODELING_META)

    if features is None or labels is None or meta is None:
        log.error("Required input files missing. Exiting.")
        sys.exit(1)

    log.info("Features: %s, Labels: %s", features.shape, labels.shape)

    # ── Part 1: Cirrhosis risk ────────────────────────────────────────────
    cirrhosis_df = cirrhosis_risk_model(features, labels, meta)
    if len(cirrhosis_df) > 0:
        out_path = os.path.join(OUT_DIR, "cirrhosis_risk_loco.csv")
        cirrhosis_df.to_csv(out_path, index=False)
        log.info("Saved: %s (%d rows)", os.path.basename(out_path), len(cirrhosis_df))
    else:
        log.warning("Cirrhosis risk model produced no output.")

    # ── Part 2: HCC molecular score ───────────────────────────────────────
    score_df, validation_df = hcc_molecular_score(features, labels, meta)
    if len(score_df) > 0:
        score_path = os.path.join(OUT_DIR, "hcc_molecular_score.csv")
        score_df.to_csv(score_path, index=False)
        log.info("Saved: %s (%d rows)", os.path.basename(score_path), len(score_df))
    if len(validation_df) > 0:
        val_path = os.path.join(OUT_DIR, "hcc_validation.csv")
        validation_df.to_csv(val_path, index=False)
        log.info("Saved: %s (%d tests)", os.path.basename(val_path), len(validation_df))

    elapsed = time.time() - t0
    log.info("")
    log.info("Done in %.1f seconds.", elapsed)


if __name__ == "__main__":
    main()
