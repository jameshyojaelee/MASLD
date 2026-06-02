#!/usr/bin/env python3
"""
151_progression_risk_score.py
Progression Risk Score (PRS) — classifies patients into S1 (stable) vs S2
(progressor) molecular subtypes and predicts continuous progression risk.

Uses LOCO-CV (leave-one-cohort-out) on F0-F2 patients with NMF subtype labels.
Five feature configurations (M1-M5) for ablation study.  Primary model is
elastic net; secondary is XGBoost (M5 only) with SHAP interpretation.

Input:
  - results/prognosis/prognosis_feature_matrix.csv  (1,444 x 253)
  - results/prognosis/prognosis_pseudo_labels.csv    (1,444 x 9)

Output (all to results/prognosis/):
  - prs_loco_results.csv           Per-fold AUROC/AUPRC per config
  - prs_all_predictions.csv        Per-sample P(S2) from LOCO-CV
  - prs_feature_importance.csv     Stable features with aggregated coefficients
  - prs_model_comparison.csv       Summary: mean AUROC +/- SD per config
  - prs_clinical_stratification.csv  Risk quartile analysis for F1-F2 patients
  - prs_external_validation.csv    Results on non-LOCO samples
  - prs_cps_regression.csv         Continuous CPS regression per fold

SLURM: io partition, 8 CPUs, 16G RAM, 48h
Env:   micromamba activate rapids_singlecell

Usage:
  sbatch --job-name=stg151_prs \
         --partition=io --cpus-per-task=8 --mem=16G --time=48:00:00 \
         --output=logs/151_prs_%j.out \
         --error=logs/151_prs_%j.err \
         --wrap="bash -c 'eval \"\\$(micromamba shell hook --shell bash)\" && \\
                 micromamba activate rapids_singlecell && \\
                 cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts && \\
                 python 151_progression_risk_score.py'"
"""

# ---------------------------------------------------------------------------
# Prevent CUDA init on CPU-only nodes
# ---------------------------------------------------------------------------
import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["RAPIDS_NO_INITIALIZE"] = "1"
os.environ["CUDF_PANDAS"] = "0"

import sys
import time
import warnings
import logging

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV, ElasticNetCV
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    brier_score_loss,
)
from xgboost import XGBClassifier

try:
    import shap
    HAS_SHAP = True
except ImportError:
    HAS_SHAP = False

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
# Suppress convergence warnings from saga solver during inner CV
warnings.filterwarnings("ignore", category=Warning, module="sklearn.linear_model")

# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------
SEED = 42
np.random.seed(SEED)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
INT = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR = os.path.join(INT, "results")
OUTDIR = os.path.join(RDIR, "prognosis")
os.makedirs(OUTDIR, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(os.path.join(OUTDIR, "151_prs.log")),
    ],
)
log = logging.getLogger(__name__)

# ===================================================================
# PART 1: Data Loading
# ===================================================================
log.info("=" * 70)
log.info("PART 1: Loading data")
log.info("=" * 70)

feat = pd.read_csv(os.path.join(OUTDIR, "prognosis_feature_matrix.csv"))
labels = pd.read_csv(os.path.join(OUTDIR, "prognosis_pseudo_labels.csv"))

# Merge on sample_id
df = labels.merge(feat, on="sample_id", how="inner")
log.info(f"Merged data: {df.shape[0]} samples x {df.shape[1]} columns")

# All feature columns (exclude label and metadata columns)
label_cols = list(labels.columns)
feature_cols = [c for c in feat.columns if c != "sample_id"]
log.info(f"Total features: {len(feature_cols)}")

# F0-F2 samples with S1/S2 labels and LOCO fold assignment
df_f02 = df[
    (df["fibrosis_stage"].isin([0, 1, 2]))
    & (df["s2_binary"].notna())
    & (df["loco_fold_fibrosis"] != "excluded")
].copy()
df_f02["s2_binary"] = df_f02["s2_binary"].astype(int)
log.info(f"F0-F2 LOCO-eligible with S2 labels: {df_f02.shape[0]} samples")
log.info(f"  S1: {(df_f02['s2_binary'] == 0).sum()}, "
         f"S2: {(df_f02['s2_binary'] == 1).sum()} "
         f"({df_f02['s2_binary'].mean()*100:.1f}% S2)")

folds = sorted(df_f02["loco_fold_fibrosis"].unique())
log.info(f"LOCO folds ({len(folds)}): {folds}")

# ===================================================================
# PART 2: Feature Configuration Ablation
# ===================================================================
log.info("")
log.info("=" * 70)
log.info("PART 2: Feature configurations")
log.info("=" * 70)

div_cols = [c for c in feature_cols if c.startswith("div_")]
ct_cols = [c for c in feature_cols if c.startswith("ct_")]
emb_cols = [c for c in feature_cols if c.startswith("nasvae_") or c.startswith("bulkformer_")]
trans_cols = [c for c in feature_cols if c.startswith("trans_")]
hcc_cols = [c for c in feature_cols if c.startswith("hcc_")]
clin_cols = [c for c in feature_cols if c.startswith("clin_")]

CONFIGS = {
    "M1_clinical": ["clin_fibrosis_stage", "clin_sex"],
    "M2_celltype": ct_cols,
    "M3_divergence": div_cols,  # will apply within-fold variance filter to top 100
    "M4_embeddings": emb_cols,
    "M5_full": feature_cols,  # all 253
}

for name, cols in CONFIGS.items():
    log.info(f"  {name}: {len(cols)} features")

# ===================================================================
# PART 3: LOCO-CV Main Loop
# ===================================================================
log.info("")
log.info("=" * 70)
log.info("PART 3: LOCO-CV Main Loop")
log.info("=" * 70)

t_start = time.time()

# Storage
loco_rows = []           # per-fold metrics
all_predictions = []     # per-sample predictions
fold_importances = {}    # config -> list of (fold, coefs) pairs
xgb_shap_values = []     # for XGBoost SHAP
cps_regression_rows = [] # continuous CPS regression


def safe_impute(X_train, X_test):
    """Impute NaN with training-set median (handles clin_age NaN)."""
    medians = np.nanmedian(X_train, axis=0)
    for j in range(X_train.shape[1]):
        if np.isnan(medians[j]):
            medians[j] = 0.0
        mask_tr = np.isnan(X_train[:, j])
        mask_te = np.isnan(X_test[:, j])
        X_train[mask_tr, j] = medians[j]
        X_test[mask_te, j] = medians[j]
    return X_train, X_test


def run_elastic_net(X_train, y_train, X_test, y_test, config_name, fold_name):
    """Run LogisticRegressionCV (elastic net) for binary S1/S2."""
    model = LogisticRegressionCV(
        penalty="elasticnet",
        l1_ratios=[0.1, 0.5, 0.9],
        solver="saga",
        cv=5,
        class_weight="balanced",
        max_iter=5000,
        scoring="roc_auc",
        random_state=SEED,
        n_jobs=-1,
    )
    model.fit(X_train, y_train)
    y_prob = model.predict_proba(X_test)[:, 1]

    auroc = roc_auc_score(y_test, y_prob) if len(np.unique(y_test)) > 1 else np.nan
    auprc = average_precision_score(y_test, y_prob) if len(np.unique(y_test)) > 1 else np.nan
    brier = brier_score_loss(y_test, y_prob)

    return model, y_prob, auroc, auprc, brier


def run_xgboost(X_train, y_train, X_test, y_test, fold_name, feat_names):
    """Run XGBoost for binary S1/S2 (M5 only)."""
    n_s1 = (y_train == 0).sum()
    n_s2 = (y_train == 1).sum()
    spw = n_s1 / max(n_s2, 1)

    model = XGBClassifier(
        n_estimators=200,
        max_depth=4,
        learning_rate=0.05,
        scale_pos_weight=spw,
        subsample=0.8,
        colsample_bytree=0.8,
        random_state=SEED,
        n_jobs=-1,
        eval_metric="logloss",
    )
    model.fit(X_train, y_train, verbose=False)
    y_prob = model.predict_proba(X_test)[:, 1]

    auroc = roc_auc_score(y_test, y_prob) if len(np.unique(y_test)) > 1 else np.nan
    auprc = average_precision_score(y_test, y_prob) if len(np.unique(y_test)) > 1 else np.nan
    brier = brier_score_loss(y_test, y_prob)

    # SHAP
    shap_vals = None
    if HAS_SHAP:
        try:
            explainer = shap.TreeExplainer(model)
            shap_vals = explainer.shap_values(X_test)
        except Exception as e:
            log.warning(f"  SHAP failed for fold {fold_name}: {e}")

    return model, y_prob, auroc, auprc, brier, shap_vals


# Main loop
for config_name, raw_cols in CONFIGS.items():
    log.info(f"\n--- Config: {config_name} ---")
    fold_importances[config_name] = []

    for fold in folds:
        train_mask = df_f02["loco_fold_fibrosis"] != fold
        test_mask = df_f02["loco_fold_fibrosis"] == fold

        train_df = df_f02[train_mask]
        test_df = df_f02[test_mask]

        if test_df.shape[0] < 5:
            log.warning(f"  Fold {fold}: too few test samples ({test_df.shape[0]}), skipping")
            continue

        y_train = train_df["s2_binary"].values
        y_test = test_df["s2_binary"].values

        # Feature selection: for M3 (divergence), use within-fold variance filter
        if config_name == "M3_divergence":
            # Select top 100 by training variance
            variances = train_df[raw_cols].var()
            top_k = min(100, len(raw_cols))
            selected_cols = variances.nlargest(top_k).index.tolist()
        else:
            selected_cols = raw_cols

        X_train = train_df[selected_cols].values.astype(np.float64)
        X_test = test_df[selected_cols].values.astype(np.float64)

        # Impute NaN
        X_train, X_test = safe_impute(X_train, X_test)

        # Scale
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        # Drop any constant columns (zero variance after scaling)
        nz_mask = X_train.std(axis=0) > 1e-10
        if not nz_mask.all():
            n_drop = (~nz_mask).sum()
            log.info(f"  Fold {fold}: dropping {n_drop} constant features")
            X_train = X_train[:, nz_mask]
            X_test = X_test[:, nz_mask]
            selected_cols = [c for c, m in zip(selected_cols, nz_mask) if m]

        log.info(f"  Fold {fold}: train {X_train.shape[0]} "
                 f"(S2={y_train.sum()}), test {X_test.shape[0]} "
                 f"(S2={y_test.sum()}), features {X_train.shape[1]}")

        # --- Elastic net ---
        en_model, en_prob, en_auroc, en_auprc, en_brier = run_elastic_net(
            X_train, y_train, X_test, y_test, config_name, fold
        )
        log.info(f"    EN: AUROC={en_auroc:.3f}  AUPRC={en_auprc:.3f}  Brier={en_brier:.3f}")

        # Store per-fold results
        loco_rows.append({
            "config": config_name,
            "fold": fold,
            "method": "elastic_net",
            "n_train": X_train.shape[0],
            "n_test": X_test.shape[0],
            "n_features": X_train.shape[1],
            "s2_prevalence_train": y_train.mean(),
            "s2_prevalence_test": y_test.mean(),
            "auroc": en_auroc,
            "auprc": en_auprc,
            "brier": en_brier,
        })

        # Store per-sample predictions
        for i, idx in enumerate(test_df.index):
            all_predictions.append({
                "sample_id": test_df.loc[idx, "sample_id"],
                "fold": fold,
                "config": config_name,
                "method": "elastic_net",
                "s2_binary_true": y_test[i],
                "p_s2_predicted": en_prob[i],
            })

        # Feature importances (elastic net coefficients)
        coefs = np.zeros(len(selected_cols))
        coefs[:] = en_model.coef_[0]
        fold_importances[config_name].append(
            pd.Series(coefs, index=selected_cols, name=fold)
        )

        # --- XGBoost (M5 only) ---
        if config_name == "M5_full":
            xgb_model, xgb_prob, xgb_auroc, xgb_auprc, xgb_brier, shap_vals = \
                run_xgboost(X_train, y_train, X_test, y_test, fold, selected_cols)
            log.info(f"    XGB: AUROC={xgb_auroc:.3f}  AUPRC={xgb_auprc:.3f}  Brier={xgb_brier:.3f}")

            loco_rows.append({
                "config": config_name,
                "fold": fold,
                "method": "xgboost",
                "n_train": X_train.shape[0],
                "n_test": X_test.shape[0],
                "n_features": X_train.shape[1],
                "s2_prevalence_train": y_train.mean(),
                "s2_prevalence_test": y_test.mean(),
                "auroc": xgb_auroc,
                "auprc": xgb_auprc,
                "brier": xgb_brier,
            })

            for i, idx in enumerate(test_df.index):
                all_predictions.append({
                    "sample_id": test_df.loc[idx, "sample_id"],
                    "fold": fold,
                    "config": config_name,
                    "method": "xgboost",
                    "s2_binary_true": y_test[i],
                    "p_s2_predicted": xgb_prob[i],
                })

            # Accumulate SHAP values
            if shap_vals is not None:
                shap_df = pd.DataFrame(shap_vals, columns=selected_cols)
                shap_df["fold"] = fold
                xgb_shap_values.append(shap_df)

        # --- Continuous CPS regression (M5 only) ---
        if config_name == "M5_full":
            cps_train = train_df["cps"].values
            cps_test = test_df["cps"].values

            # Check for NaN in CPS
            valid_tr = ~np.isnan(cps_train)
            valid_te = ~np.isnan(cps_test)

            if valid_tr.sum() >= 10 and valid_te.sum() >= 5:
                reg = ElasticNetCV(
                    l1_ratio=[0.1, 0.5, 0.9],
                    cv=5,
                    random_state=SEED,
                    max_iter=5000,
                    n_jobs=-1,
                )
                reg.fit(X_train[valid_tr], cps_train[valid_tr])
                cps_pred = reg.predict(X_test[valid_te])
                rho, pval = spearmanr(cps_test[valid_te], cps_pred)
                mae = np.mean(np.abs(cps_test[valid_te] - cps_pred))
                log.info(f"    CPS regression: rho={rho:.3f} (p={pval:.2e}), MAE={mae:.3f}")

                cps_regression_rows.append({
                    "fold": fold,
                    "n_train": valid_tr.sum(),
                    "n_test": valid_te.sum(),
                    "spearman_rho": rho,
                    "spearman_pval": pval,
                    "mae": mae,
                    "l1_ratio_best": reg.l1_ratio_,
                    "alpha_best": reg.alpha_,
                })

elapsed = time.time() - t_start
log.info(f"\nLOCO-CV complete in {elapsed/60:.1f} minutes")

# ===================================================================
# PART 4: External Validation
# ===================================================================
log.info("")
log.info("=" * 70)
log.info("PART 4: External Validation")
log.info("=" * 70)

# Train M5 on all LOCO-eligible F0-F2 samples, predict on excluded samples
# Excluded = samples with s2_binary but no LOCO fold (GSE213621 etc.) OR
# samples with fibrosis staging but excluded from LOCO
ext_train = df_f02.copy()  # 492 LOCO-eligible F0-F2

# External: samples with s2_binary labels but NOT in LOCO F0-F2 set
ext_test = df[
    (df["s2_binary"].notna())
    & (~df.index.isin(ext_train.index))
].copy()
ext_test["s2_binary"] = ext_test["s2_binary"].astype(int)

log.info(f"External train: {ext_train.shape[0]} samples (LOCO-eligible F0-F2)")
log.info(f"External test: {ext_test.shape[0]} samples (non-LOCO with S2 labels)")

ext_results = []
if ext_test.shape[0] >= 10:
    X_ext_train = ext_train[feature_cols].values.astype(np.float64)
    X_ext_test = ext_test[feature_cols].values.astype(np.float64)
    y_ext_train = ext_train["s2_binary"].values
    y_ext_test = ext_test["s2_binary"].values

    X_ext_train, X_ext_test = safe_impute(X_ext_train, X_ext_test)

    scaler_ext = StandardScaler()
    X_ext_train = scaler_ext.fit_transform(X_ext_train)
    X_ext_test = scaler_ext.transform(X_ext_test)

    # Drop constant
    nz = X_ext_train.std(axis=0) > 1e-10
    X_ext_train = X_ext_train[:, nz]
    X_ext_test = X_ext_test[:, nz]

    # Elastic net
    ext_en = LogisticRegressionCV(
        penalty="elasticnet",
        l1_ratios=[0.1, 0.5, 0.9],
        solver="saga",
        cv=5,
        class_weight="balanced",
        max_iter=5000,
        scoring="roc_auc",
        random_state=SEED,
        n_jobs=-1,
    )
    ext_en.fit(X_ext_train, y_ext_train)
    ext_prob = ext_en.predict_proba(X_ext_test)[:, 1]

    ext_auroc = roc_auc_score(y_ext_test, ext_prob) if len(np.unique(y_ext_test)) > 1 else np.nan
    ext_auprc = average_precision_score(y_ext_test, ext_prob) if len(np.unique(y_ext_test)) > 1 else np.nan
    ext_brier = brier_score_loss(y_ext_test, ext_prob)

    log.info(f"  External EN: AUROC={ext_auroc:.3f}  AUPRC={ext_auprc:.3f}  Brier={ext_brier:.3f}")
    log.info(f"  S2 prevalence: train={y_ext_train.mean():.3f}, test={y_ext_test.mean():.3f}")

    # Per-dataset breakdown
    for ds in sorted(ext_test["dataset"].unique()):
        ds_mask = ext_test["dataset"] == ds
        ds_y = y_ext_test[ds_mask.values]
        ds_p = ext_prob[ds_mask.values]
        ds_auroc = roc_auc_score(ds_y, ds_p) if len(np.unique(ds_y)) > 1 else np.nan
        ds_auprc = average_precision_score(ds_y, ds_p) if len(np.unique(ds_y)) > 1 else np.nan
        log.info(f"    {ds}: n={ds_mask.sum()}, S2={ds_y.sum()}, AUROC={ds_auroc:.3f}")

        ext_results.append({
            "dataset": ds,
            "n_samples": int(ds_mask.sum()),
            "n_s2": int(ds_y.sum()),
            "s2_prevalence": ds_y.mean(),
            "auroc": ds_auroc,
            "auprc": ds_auprc,
        })

    # Overall
    ext_results.append({
        "dataset": "ALL_EXTERNAL",
        "n_samples": ext_test.shape[0],
        "n_s2": int(y_ext_test.sum()),
        "s2_prevalence": y_ext_test.mean(),
        "auroc": ext_auroc,
        "auprc": ext_auprc,
    })
else:
    log.info("  Insufficient external samples; skipping.")

# ===================================================================
# PART 5: Clinical Stratification
# ===================================================================
log.info("")
log.info("=" * 70)
log.info("PART 5: Clinical Stratification (F1-F2 patients)")
log.info("=" * 70)

# Get M5 elastic net predictions from LOCO-CV
preds_df = pd.DataFrame(all_predictions)
m5_en_preds = preds_df[
    (preds_df["config"] == "M5_full") & (preds_df["method"] == "elastic_net")
].copy()

# Merge back to get clinical info
m5_en_preds = m5_en_preds.merge(
    df[["sample_id", "fibrosis_stage", "nas_score", "dataset"]],
    on="sample_id",
    how="left",
)

# Also merge cell-type fractions for quartile characterization
ct_feat_cols = [c for c in feature_cols if c.startswith("ct_")]
m5_en_preds = m5_en_preds.merge(
    feat[["sample_id"] + ct_feat_cols],
    on="sample_id",
    how="left",
)

# Focus on F1-F2
f12 = m5_en_preds[m5_en_preds["fibrosis_stage"].isin([1, 2])].copy()
log.info(f"F1-F2 patients with LOCO predictions: {f12.shape[0]}")

strat_rows = []
if f12.shape[0] >= 20:
    # Risk quartiles
    f12["risk_quartile"] = pd.qcut(
        f12["p_s2_predicted"], q=4, labels=["Q1_low", "Q2", "Q3", "Q4_high"]
    )

    for q in ["Q1_low", "Q2", "Q3", "Q4_high"]:
        sub = f12[f12["risk_quartile"] == q]
        row = {
            "quartile": q,
            "n_samples": len(sub),
            "n_s2_actual": int(sub["s2_binary_true"].sum()),
            "s2_fraction": sub["s2_binary_true"].mean(),
            "mean_p_s2": sub["p_s2_predicted"].mean(),
            "median_fibrosis": sub["fibrosis_stage"].median(),
            "median_stellate": sub["ct_Stellate"].median() if "ct_Stellate" in sub.columns else np.nan,
            "mean_immune_infiltration": sub["ct_immune_infiltration"].mean() if "ct_immune_infiltration" in sub.columns else np.nan,
        }
        if sub["nas_score"].notna().sum() > 0:
            row["median_nas"] = sub["nas_score"].median()
        else:
            row["median_nas"] = np.nan
        strat_rows.append(row)
        log.info(f"  {q}: n={row['n_samples']}, S2={row['n_s2_actual']} "
                 f"({row['s2_fraction']*100:.1f}%), "
                 f"median_fib={row['median_fibrosis']}, "
                 f"median_stellate={row['median_stellate']:.4f}")

    # Net reclassification improvement (NRI)
    # Baseline: fibrosis stage alone (F2 = high risk, F1 = low risk)
    f12["baseline_risk"] = (f12["fibrosis_stage"] >= 2).astype(int)
    f12["prs_risk"] = (f12["p_s2_predicted"] >= 0.5).astype(int)
    y_true = f12["s2_binary_true"].values

    # Events (S2=1): fraction reclassified up by PRS
    events = f12[f12["s2_binary_true"] == 1]
    nonevents = f12[f12["s2_binary_true"] == 0]

    if len(events) > 0 and len(nonevents) > 0:
        # NRI for events: gained correct classification
        event_up = ((events["prs_risk"] == 1) & (events["baseline_risk"] == 0)).sum()
        event_down = ((events["prs_risk"] == 0) & (events["baseline_risk"] == 1)).sum()
        nri_events = (event_up - event_down) / len(events)

        nonevent_down = ((nonevents["prs_risk"] == 0) & (nonevents["baseline_risk"] == 1)).sum()
        nonevent_up = ((nonevents["prs_risk"] == 1) & (nonevents["baseline_risk"] == 0)).sum()
        nri_nonevents = (nonevent_down - nonevent_up) / len(nonevents)

        nri_total = nri_events + nri_nonevents
        log.info(f"  NRI (events): {nri_events:.3f}")
        log.info(f"  NRI (non-events): {nri_nonevents:.3f}")
        log.info(f"  NRI (total): {nri_total:.3f}")

        strat_rows.append({
            "quartile": "NRI_SUMMARY",
            "n_samples": len(f12),
            "n_s2_actual": int(y_true.sum()),
            "s2_fraction": y_true.mean(),
            "mean_p_s2": nri_total,
            "median_fibrosis": nri_events,
            "median_stellate": nri_nonevents,
            "mean_immune_infiltration": np.nan,
            "median_nas": np.nan,
        })

    # Decision curve: net benefit at thresholds
    log.info("  Decision curve analysis:")
    for thresh in [0.1, 0.2, 0.3, 0.4, 0.5]:
        pred_pos = (f12["p_s2_predicted"] >= thresh).values
        tp = ((pred_pos) & (y_true == 1)).sum()
        fp = ((pred_pos) & (y_true == 0)).sum()
        n = len(y_true)
        net_benefit = tp / n - fp / n * (thresh / (1 - thresh)) if thresh < 1 else 0
        log.info(f"    threshold={thresh:.1f}: net_benefit={net_benefit:.4f}")

        strat_rows.append({
            "quartile": f"DC_thresh_{thresh}",
            "n_samples": n,
            "n_s2_actual": int(y_true.sum()),
            "s2_fraction": thresh,
            "mean_p_s2": net_benefit,
            "median_fibrosis": int(tp),
            "median_stellate": int(fp),
            "mean_immune_infiltration": np.nan,
            "median_nas": np.nan,
        })
else:
    log.warning("  Insufficient F1-F2 patients for stratification analysis")

# ===================================================================
# PART 6: Feature Importance
# ===================================================================
log.info("")
log.info("=" * 70)
log.info("PART 6: Feature Importance")
log.info("=" * 70)

importance_rows = []

for config_name, coef_list in fold_importances.items():
    if not coef_list:
        continue

    # Aggregate coefficients: align on feature names
    all_feats = set()
    for s in coef_list:
        all_feats.update(s.index)
    all_feats = sorted(all_feats)

    coef_matrix = pd.DataFrame(index=all_feats, columns=[s.name for s in coef_list])
    for s in coef_list:
        coef_matrix.loc[s.index, s.name] = s.values

    coef_matrix = coef_matrix.astype(float)

    # Mean coefficient and stability (non-zero in how many folds)
    mean_coef = coef_matrix.mean(axis=1)
    nonzero_count = (coef_matrix.abs() > 1e-10).sum(axis=1)
    stability = nonzero_count / len(coef_list)

    for feat_name in all_feats:
        importance_rows.append({
            "config": config_name,
            "feature": feat_name,
            "mean_coefficient": mean_coef[feat_name],
            "abs_mean_coefficient": abs(mean_coef[feat_name]),
            "n_folds_nonzero": int(nonzero_count[feat_name]),
            "n_folds_total": len(coef_list),
            "stability": stability[feat_name],
        })

# SHAP for XGBoost M5
if xgb_shap_values:
    log.info("  Aggregating XGBoost SHAP values...")
    shap_all = pd.concat(xgb_shap_values, ignore_index=True)
    shap_cols = [c for c in shap_all.columns if c != "fold"]
    mean_abs_shap = shap_all[shap_cols].abs().mean()
    top_shap = mean_abs_shap.nlargest(20)

    log.info("  Top 20 XGBoost features (mean |SHAP|):")
    for feat_name, val in top_shap.items():
        log.info(f"    {feat_name}: {val:.4f}")
        importance_rows.append({
            "config": "M5_full_xgboost_shap",
            "feature": feat_name,
            "mean_coefficient": val,  # store SHAP as "coefficient" for unified format
            "abs_mean_coefficient": val,
            "n_folds_nonzero": len(xgb_shap_values),
            "n_folds_total": len(xgb_shap_values),
            "stability": 1.0,
        })

# Identify stable features for M5 elastic net
m5_imp = [r for r in importance_rows
          if r["config"] == "M5_full" and r["stability"] >= 4/6]
m5_imp_sorted = sorted(m5_imp, key=lambda x: x["abs_mean_coefficient"], reverse=True)

log.info(f"\nStable M5 EN features (non-zero in >=4/6 folds): {len(m5_imp_sorted)}")
log.info("Top 20 stable features:")
for r in m5_imp_sorted[:20]:
    log.info(f"  {r['feature']}: coef={r['mean_coefficient']:.4f}, "
             f"stability={r['n_folds_nonzero']}/{r['n_folds_total']}")

# ===================================================================
# PART 7: Save Results
# ===================================================================
log.info("")
log.info("=" * 70)
log.info("PART 7: Saving results")
log.info("=" * 70)

# 1. Per-fold results
loco_df = pd.DataFrame(loco_rows)
loco_df.to_csv(os.path.join(OUTDIR, "prs_loco_results.csv"), index=False)
log.info(f"  Saved prs_loco_results.csv ({loco_df.shape[0]} rows)")

# 2. All predictions
pred_df = pd.DataFrame(all_predictions)
pred_df.to_csv(os.path.join(OUTDIR, "prs_all_predictions.csv"), index=False)
log.info(f"  Saved prs_all_predictions.csv ({pred_df.shape[0]} rows)")

# 3. Feature importance
imp_df = pd.DataFrame(importance_rows)
imp_df.to_csv(os.path.join(OUTDIR, "prs_feature_importance.csv"), index=False)
log.info(f"  Saved prs_feature_importance.csv ({imp_df.shape[0]} rows)")

# 4. Model comparison summary
summary_rows = []
for config_name in CONFIGS:
    for method in ["elastic_net", "xgboost"]:
        sub = loco_df[(loco_df["config"] == config_name) & (loco_df["method"] == method)]
        if sub.shape[0] == 0:
            continue
        summary_rows.append({
            "config": config_name,
            "method": method,
            "n_folds": sub.shape[0],
            "mean_auroc": sub["auroc"].mean(),
            "sd_auroc": sub["auroc"].std(),
            "mean_auprc": sub["auprc"].mean(),
            "sd_auprc": sub["auprc"].std(),
            "mean_brier": sub["brier"].mean(),
            "sd_brier": sub["brier"].std(),
            "min_auroc": sub["auroc"].min(),
            "max_auroc": sub["auroc"].max(),
        })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(OUTDIR, "prs_model_comparison.csv"), index=False)
log.info(f"  Saved prs_model_comparison.csv ({summary_df.shape[0]} rows)")

# Print summary table
log.info("\n=== MODEL COMPARISON SUMMARY ===")
for _, row in summary_df.iterrows():
    log.info(f"  {row['config']:20s} {row['method']:12s}  "
             f"AUROC={row['mean_auroc']:.3f}+/-{row['sd_auroc']:.3f}  "
             f"AUPRC={row['mean_auprc']:.3f}+/-{row['sd_auprc']:.3f}  "
             f"Brier={row['mean_brier']:.3f}")

# 5. Clinical stratification
if strat_rows:
    strat_df = pd.DataFrame(strat_rows)
    strat_df.to_csv(os.path.join(OUTDIR, "prs_clinical_stratification.csv"), index=False)
    log.info(f"  Saved prs_clinical_stratification.csv ({strat_df.shape[0]} rows)")

# 6. External validation
if ext_results:
    ext_df = pd.DataFrame(ext_results)
    ext_df.to_csv(os.path.join(OUTDIR, "prs_external_validation.csv"), index=False)
    log.info(f"  Saved prs_external_validation.csv ({ext_df.shape[0]} rows)")

# 7. CPS regression
if cps_regression_rows:
    cps_df = pd.DataFrame(cps_regression_rows)
    cps_df.to_csv(os.path.join(OUTDIR, "prs_cps_regression.csv"), index=False)
    log.info(f"  Saved prs_cps_regression.csv ({cps_df.shape[0]} rows)")
    log.info(f"  CPS mean Spearman rho: {cps_df['spearman_rho'].mean():.3f} "
             f"+/- {cps_df['spearman_rho'].std():.3f}")

# Final timing
total_time = time.time() - t_start
log.info(f"\nTotal runtime: {total_time/60:.1f} minutes")
log.info("Script 151 complete.")
