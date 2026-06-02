#!/usr/bin/env python3
"""
108_stacking_ensemble.py
Stacking meta-learner that combines predictions from all base models.

Base models (load out-of-fold predictions from existing result files):
  1. V3 elastic net          (v3_proper_cv_results.csv, target=fib_ge3)
  2. NAS-VAE + RandomForest  (nas_vae_classical_results.csv)
  3. Concept bottleneck       (concept_bottleneck_results.csv)
  4. Multi-task CORN          (multitask_nas_results.csv)

Stacking targets:
  - F>=3 binary: stack 4 model probabilities -> logistic regression meta-learner
    (leave-one-fold-out)
  - NAS ordinal: stack available model QWK/class predictions -> ordinal meta-learner

Output (all to results/staging_classifier/):
  - stacking_ensemble_results.csv         (per-sample stacking predictions)
  - stacking_vs_single_comparison.csv     (stacking AUROC vs best single model)

SLURM: cpu partition, 8 CPUs, 32G RAM, 48h
Env:   micromamba activate rapids_singlecell

Usage:
  sbatch --job-name=stg108_stacking \
         --partition=cpu --cpus-per-task=8 --mem=32G --time=48:00:00 \
         --output=logs/108_stacking_%j.out \
         --error=logs/108_stacking_%j.err \
         --wrap="bash -c 'eval \"\\$(micromamba shell hook --shell bash)\" && \\
                 micromamba activate rapids_singlecell && \\
                 cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts && \\
                 python 108_stacking_ensemble.py'"
"""

# Must be set before ANY other imports to prevent CUDA crash on cpu nodes
import os
os.environ["RAPIDS_NO_INITIALIZE"] = "1"
os.environ["CUDF_PANDAS"] = "0"

import sys
import time
import warnings
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, accuracy_score, cohen_kappa_score
from scipy import stats

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

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
OUTDIR = os.path.join(RDIR, "staging_classifier")
LOGDIR = os.path.join(OUTDIR, "logs")
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(LOGDIR, exist_ok=True)

# Input files (out-of-fold predictions from base models)
V3_RESULTS = os.path.join(OUTDIR, "v3_proper_cv_results.csv")
VAE_RESULTS = os.path.join(OUTDIR, "nas_vae_classical_results.csv")
CBM_RESULTS = os.path.join(OUTDIR, "concept_bottleneck_results.csv")
MTL_RESULTS = os.path.join(OUTDIR, "multitask_nas_results.csv")
META_PATH = os.path.join(OUTDIR, "modeling_metadata.csv")

# Summary files for single-model comparison
V3_SUMMARY = os.path.join(OUTDIR, "v3_model_summary.csv")
VAE_SUMMARY = os.path.join(OUTDIR, "nas_vae_model_summary.csv")
CBM_SUMMARY = os.path.join(OUTDIR, "concept_bottleneck_summary.csv")
MTL_SUMMARY = os.path.join(OUTDIR, "multitask_nas_summary.csv")

# Outputs
OUT_RESULTS = os.path.join(OUTDIR, "stacking_ensemble_results.csv")
OUT_COMPARISON = os.path.join(OUTDIR, "stacking_vs_single_comparison.csv")

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(os.path.join(LOGDIR, "108_stacking_ensemble.log")),
    ],
)
log = logging.getLogger(__name__)

print("=" * 70)
print("108: Stacking Ensemble Meta-Learner for MASLD Staging")
print("=" * 70)


# ===================================================================
# 1. Load out-of-fold predictions from base models
# ===================================================================

def load_v3_predictions():
    """
    Load per-sample out-of-fold predictions from V3 elastic net.
    Format: sample_id, dataset, true_label, predicted, probability, fold, target
    """
    if not os.path.exists(V3_RESULTS):
        log.warning("V3 results not found: %s", V3_RESULTS)
        return None

    df = pd.read_csv(V3_RESULTS)
    log.info("V3 predictions loaded: %d rows, targets: %s",
             len(df), df["target"].unique().tolist())
    return df


def load_vae_classical_predictions():
    """
    Load NAS-VAE + classical ML results.
    Format: target, target_desc, fold, model, n_train, n_test, n_classes,
            accuracy, balanced_accuracy, auroc, f1_macro, kappa, qwk
    These are fold-level summaries, not per-sample. We extract fold-level
    AUROC for the best model (RandomForest) per target.
    """
    if not os.path.exists(VAE_RESULTS):
        log.warning("VAE classical results not found: %s", VAE_RESULTS)
        return None

    df = pd.read_csv(VAE_RESULTS)
    log.info("VAE classical results loaded: %d rows, targets: %s",
             len(df), df["target"].unique().tolist())
    return df


def load_cbm_predictions():
    """
    Load concept bottleneck model results.
    Format: task, n_valid, auroc, acc, fold, fold_type, qwk, mae, adj_acc
    """
    if not os.path.exists(CBM_RESULTS):
        log.warning("Concept bottleneck results not found: %s", CBM_RESULTS)
        return None

    df = pd.read_csv(CBM_RESULTS)
    log.info("CBM results loaded: %d rows, tasks: %s",
             len(df), df["task"].unique().tolist())
    return df


def load_mtl_predictions():
    """
    Load multi-task CORN results.
    Format: task, n_valid, qwk, mae, acc, adj_acc, fold, auroc
    """
    if not os.path.exists(MTL_RESULTS):
        log.warning("Multi-task results not found: %s", MTL_RESULTS)
        return None

    df = pd.read_csv(MTL_RESULTS)
    log.info("MTL results loaded: %d rows, tasks: %s",
             len(df), df["task"].unique().tolist())
    return df


def load_metadata():
    """Load modeling metadata for ground-truth labels."""
    meta = pd.read_csv(META_PATH)
    log.info("Metadata loaded: %d samples, columns: %s",
             len(meta), meta.columns.tolist()[:10])
    return meta


# ===================================================================
# 2. Build per-sample stacking matrix for F>=3 binary
# ===================================================================

def build_fib_ge3_stacking_matrix(v3_df, meta):
    """
    Build per-sample stacking matrix for F>=3 binary classification.

    V3 elastic net provides per-sample out-of-fold probabilities.
    For models without per-sample predictions, we use leave-one-cohort-out
    (LOCO) fold-level embeddings to reconstruct per-sample features.

    Returns: DataFrame with columns [sample_id, fold, y_true, p_v3, p_vae, p_cbm, p_mtl]
    """
    # --- V3: direct per-sample probabilities ---
    v3_fib = v3_df[v3_df["target"] == "fib_ge3"].copy()
    if v3_fib.empty:
        log.error("No fib_ge3 target in V3 results")
        return None

    # Build sample-level matrix
    stack = v3_fib[["sample_id", "fold", "true_label", "probability"]].rename(
        columns={"true_label": "y_true", "probability": "p_v3", "fold": "loco_fold"}
    )

    log.info("Stacking matrix for fib_ge3: %d samples from V3", len(stack))

    # --- VAE classical: reconstruct per-sample probabilities from LOCO embeddings ---
    # Load LOCO NAS embeddings + retrain RF for each fold
    p_vae = _reconstruct_vae_per_sample_probs(meta, target="fib_ge3")
    if p_vae is not None:
        stack = stack.merge(p_vae, on="sample_id", how="left")
    else:
        stack["p_vae"] = np.nan

    # --- CBM: reconstruct from concept activations + metadata ---
    p_cbm = _reconstruct_cbm_per_sample_probs(meta, target="fib_ge3")
    if p_cbm is not None:
        stack = stack.merge(p_cbm, on="sample_id", how="left")
    else:
        stack["p_cbm"] = np.nan

    # --- MTL: reconstruct from LOCO embeddings + NAS VAE ---
    p_mtl = _reconstruct_mtl_per_sample_probs(meta, target="fib_ge3")
    if p_mtl is not None:
        stack = stack.merge(p_mtl, on="sample_id", how="left")
    else:
        stack["p_mtl"] = np.nan

    # Keep only samples with valid labels
    stack = stack[stack["y_true"] >= 0].copy()
    log.info("Stacking matrix after label filter: %d samples", len(stack))

    return stack


def _reconstruct_vae_per_sample_probs(meta, target="fib_ge3"):
    """
    Reconstruct per-sample probabilities from NAS-VAE embeddings using
    leave-one-cohort-out cross-validation with RandomForest.
    """
    from sklearn.ensemble import RandomForestClassifier

    embed_path = os.path.join(OUTDIR, "nas_embeddings_all_samples.csv")
    if not os.path.exists(embed_path):
        log.warning("NAS embeddings not found for VAE reconstruction")
        return None

    embeddings = pd.read_csv(embed_path)
    z_cols = [c for c in embeddings.columns if c.startswith("z")]

    # Merge with metadata for labels
    merged = embeddings.merge(meta[["sample_id", "dataset", target]], on="sample_id")
    merged = merged[merged[target] >= 0].copy()

    if len(merged) == 0:
        return None

    folds = merged["dataset"].unique()
    results = []

    for fold_ds in folds:
        train_mask = merged["dataset"] != fold_ds
        test_mask = merged["dataset"] == fold_ds

        X_train = merged.loc[train_mask, z_cols].values
        y_train = merged.loc[train_mask, target].values
        X_test = merged.loc[test_mask, z_cols].values
        test_ids = merged.loc[test_mask, "sample_id"].values

        if len(np.unique(y_train)) < 2 or test_mask.sum() == 0:
            continue

        rf = RandomForestClassifier(
            n_estimators=500, max_depth=10, min_samples_leaf=5,
            class_weight="balanced", random_state=SEED, n_jobs=-1
        )
        rf.fit(X_train, y_train)
        probs = rf.predict_proba(X_test)
        # Probability of positive class
        p_pos = probs[:, 1] if probs.shape[1] == 2 else probs[:, -1]

        for sid, p in zip(test_ids, p_pos):
            results.append({"sample_id": sid, "p_vae": float(p)})

    if not results:
        return None

    return pd.DataFrame(results)


def _reconstruct_cbm_per_sample_probs(meta, target="fib_ge3"):
    """
    Reconstruct per-sample probabilities from concept activations using
    leave-one-cohort-out logistic regression on concept features.
    """
    concept_path = os.path.join(OUTDIR, "concept_activations.csv")
    if not os.path.exists(concept_path):
        log.warning("Concept activations not found for CBM reconstruction")
        return None

    concepts = pd.read_csv(concept_path)
    c_cols = [c for c in concepts.columns if c != "sample_id"]

    merged = concepts.merge(meta[["sample_id", "dataset", target]], on="sample_id")
    merged = merged[merged[target] >= 0].copy()

    if len(merged) == 0:
        return None

    folds = merged["dataset"].unique()
    results = []

    for fold_ds in folds:
        train_mask = merged["dataset"] != fold_ds
        test_mask = merged["dataset"] == fold_ds

        X_train = merged.loc[train_mask, c_cols].values
        y_train = merged.loc[train_mask, target].values
        X_test = merged.loc[test_mask, c_cols].values
        test_ids = merged.loc[test_mask, "sample_id"].values

        if len(np.unique(y_train)) < 2 or test_mask.sum() == 0:
            continue

        lr = LogisticRegression(
            penalty="l2", C=1.0, max_iter=2000, solver="lbfgs",
            class_weight="balanced", random_state=SEED
        )
        lr.fit(X_train, y_train)
        probs = lr.predict_proba(X_test)
        p_pos = probs[:, 1] if probs.shape[1] == 2 else probs[:, -1]

        for sid, p in zip(test_ids, p_pos):
            results.append({"sample_id": sid, "p_cbm": float(p)})

    if not results:
        return None

    return pd.DataFrame(results)


def _reconstruct_mtl_per_sample_probs(meta, target="fib_ge3"):
    """
    Reconstruct per-sample probabilities for the multi-task model using
    NAS embeddings with XGBoost (mirroring CORN's feature representation).
    """
    try:
        import xgboost as xgb
    except ImportError:
        log.warning("XGBoost not available; skipping MTL reconstruction")
        return None

    embed_path = os.path.join(OUTDIR, "nas_embeddings_all_samples.csv")
    if not os.path.exists(embed_path):
        log.warning("NAS embeddings not found for MTL reconstruction")
        return None

    embeddings = pd.read_csv(embed_path)
    z_cols = [c for c in embeddings.columns if c.startswith("z")]

    # Also incorporate concept activations if available
    concept_path = os.path.join(OUTDIR, "concept_activations.csv")
    if os.path.exists(concept_path):
        concepts = pd.read_csv(concept_path)
        embeddings = embeddings.merge(concepts, on="sample_id", how="left")
        c_cols = [c for c in concepts.columns if c != "sample_id"]
        feat_cols = z_cols + c_cols
    else:
        feat_cols = z_cols

    merged = embeddings.merge(meta[["sample_id", "dataset", target]], on="sample_id")
    merged = merged[merged[target] >= 0].copy()

    if len(merged) == 0:
        return None

    folds = merged["dataset"].unique()
    results = []

    for fold_ds in folds:
        train_mask = merged["dataset"] != fold_ds
        test_mask = merged["dataset"] == fold_ds

        X_train = merged.loc[train_mask, feat_cols].values
        y_train = merged.loc[train_mask, target].values
        X_test = merged.loc[test_mask, feat_cols].values
        test_ids = merged.loc[test_mask, "sample_id"].values

        if len(np.unique(y_train)) < 2 or test_mask.sum() == 0:
            continue

        model = xgb.XGBClassifier(
            n_estimators=300, max_depth=6, learning_rate=0.05,
            subsample=0.8, colsample_bytree=0.8,
            scale_pos_weight=float(np.sum(y_train == 0)) / max(np.sum(y_train == 1), 1),
            random_state=SEED, n_jobs=-1,
            eval_metric="logloss"
        )
        model.fit(X_train, y_train, verbose=False)
        probs = model.predict_proba(X_test)
        p_pos = probs[:, 1] if probs.shape[1] == 2 else probs[:, -1]

        for sid, p in zip(test_ids, p_pos):
            results.append({"sample_id": sid, "p_mtl": float(p)})

    if not results:
        return None

    return pd.DataFrame(results)


# ===================================================================
# 3. Stacking meta-learner (leave-one-fold-out)
# ===================================================================

def run_fib_ge3_stacking(stack_df):
    """
    Train stacking meta-learner for F>=3 binary using leave-one-fold-out.

    Meta-learner: Logistic regression on base model probabilities.
    """
    prob_cols = [c for c in ["p_v3", "p_vae", "p_cbm", "p_mtl"]
                 if c in stack_df.columns]

    # Drop columns that are all NaN
    valid_cols = [c for c in prob_cols if stack_df[c].notna().sum() > 0]
    if not valid_cols:
        log.error("No valid probability columns for stacking")
        return None

    log.info("Stacking with %d base models: %s", len(valid_cols), valid_cols)

    # Fill remaining NaN with 0.5 (neutral probability)
    for c in valid_cols:
        stack_df[c] = stack_df[c].fillna(0.5)

    folds = stack_df["loco_fold"].unique()
    all_preds = []

    for fold in folds:
        train_mask = stack_df["loco_fold"] != fold
        test_mask = stack_df["loco_fold"] == fold

        X_train = stack_df.loc[train_mask, valid_cols].values
        y_train = stack_df.loc[train_mask, "y_true"].values
        X_test = stack_df.loc[test_mask, valid_cols].values
        test_ids = stack_df.loc[test_mask, "sample_id"].values
        y_test = stack_df.loc[test_mask, "y_true"].values

        if len(np.unique(y_train)) < 2 or test_mask.sum() == 0:
            log.warning("Fold %s skipped: <2 classes in train or empty test", fold)
            continue

        meta_lr = LogisticRegression(
            penalty="l2", C=1.0, max_iter=2000, solver="lbfgs",
            class_weight="balanced", random_state=SEED
        )
        meta_lr.fit(X_train, y_train)

        probs = meta_lr.predict_proba(X_test)
        p_stack = probs[:, 1] if probs.shape[1] == 2 else probs[:, -1]
        preds = meta_lr.predict(X_test)

        for sid, yt, pred, ps in zip(test_ids, y_test, preds, p_stack):
            all_preds.append({
                "sample_id": sid,
                "fold": fold,
                "y_true": int(yt),
                "y_pred_stack": int(pred),
                "p_stack": float(ps),
            })

        # Per-fold AUROC
        if len(np.unique(y_test)) >= 2:
            fold_auroc = roc_auc_score(y_test, p_stack)
            log.info("  Fold %s: AUROC = %.4f (n=%d)", fold, fold_auroc, len(y_test))

    if not all_preds:
        return None

    return pd.DataFrame(all_preds)


# ===================================================================
# 4. NAS ordinal stacking
# ===================================================================

def build_nas_ordinal_stacking_matrix(meta):
    """
    Build stacking matrix for NAS ordinal prediction.
    Uses NAS-VAE embeddings + concept activations as shared base features.
    Trains multiple ordinal classifiers and stacks their QWK predictions.
    """
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier

    embed_path = os.path.join(OUTDIR, "nas_embeddings_all_samples.csv")
    concept_path = os.path.join(OUTDIR, "concept_activations.csv")

    if not os.path.exists(embed_path):
        log.warning("NAS embeddings not found for ordinal stacking")
        return None

    embeddings = pd.read_csv(embed_path)
    z_cols = [c for c in embeddings.columns if c.startswith("z")]
    feat_cols = list(z_cols)

    if os.path.exists(concept_path):
        concepts = pd.read_csv(concept_path)
        embeddings = embeddings.merge(concepts, on="sample_id", how="left")
        c_cols = [c for c in concepts.columns if c != "sample_id"]
        feat_cols += c_cols

    # Merge with metadata
    target = "nas_score"
    merged = embeddings.merge(
        meta[["sample_id", "dataset", "nas_score", "nas_group4"]],
        on="sample_id"
    )

    # NAS 4-group ordinal: 0-1, 2-3, 4-5, 6-8
    merged = merged[merged["nas_group4"] >= 0].copy()
    if len(merged) == 0:
        log.warning("No valid NAS samples for ordinal stacking")
        return None

    log.info("NAS ordinal stacking: %d samples", len(merged))

    # LOCO folds based on datasets that have NAS scores
    folds = merged["dataset"].unique()
    all_preds = []

    for fold_ds in folds:
        train_mask = merged["dataset"] != fold_ds
        test_mask = merged["dataset"] == fold_ds

        X_train = merged.loc[train_mask, feat_cols].values
        y_train = merged.loc[train_mask, "nas_group4"].values.astype(int)
        X_test = merged.loc[test_mask, feat_cols].values
        y_test = merged.loc[test_mask, "nas_group4"].values.astype(int)
        test_ids = merged.loc[test_mask, "sample_id"].values

        if len(np.unique(y_train)) < 2 or test_mask.sum() == 0:
            continue

        # Train 3 base learners
        models = {
            "RF": RandomForestClassifier(
                n_estimators=500, max_depth=10, min_samples_leaf=5,
                class_weight="balanced", random_state=SEED, n_jobs=-1
            ),
            "GB": GradientBoostingClassifier(
                n_estimators=300, max_depth=5, learning_rate=0.05,
                subsample=0.8, random_state=SEED
            ),
            "LR": LogisticRegression(
                penalty="l2", C=1.0, max_iter=2000, solver="lbfgs",
                class_weight="balanced", random_state=SEED
            ),
        }

        # Inner LOCO for stacking (but since we only have 1 held-out fold,
        # use the OOF predictions directly from inner 5-fold CV on train set)
        from sklearn.model_selection import StratifiedKFold

        inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
        oof_train = np.zeros((X_train.shape[0], len(models)))
        test_probs = np.zeros((X_test.shape[0], len(models)))

        for m_idx, (m_name, model) in enumerate(models.items()):
            # Inner OOF for train
            inner_oof = np.zeros(X_train.shape[0])
            inner_test_preds = np.zeros(X_test.shape[0])
            n_inner = 0

            for inner_train_idx, inner_val_idx in inner_cv.split(X_train, y_train):
                Xi_train = X_train[inner_train_idx]
                yi_train = y_train[inner_train_idx]
                Xi_val = X_train[inner_val_idx]

                model_clone = type(model)(**model.get_params())
                model_clone.fit(Xi_train, yi_train)

                # Use predicted class (ordinal) as stacking feature
                inner_oof[inner_val_idx] = model_clone.predict(Xi_val).astype(float)
                inner_test_preds += model_clone.predict(X_test).astype(float)
                n_inner += 1

            oof_train[:, m_idx] = inner_oof
            test_probs[:, m_idx] = inner_test_preds / max(n_inner, 1)

        # Meta-learner: logistic regression on stacking features
        meta_lr = LogisticRegression(
            penalty="l2", C=1.0, max_iter=2000, solver="lbfgs",
            class_weight="balanced", random_state=SEED
        )
        meta_lr.fit(oof_train, y_train)
        meta_preds = meta_lr.predict(test_probs)

        # QWK for this fold
        if len(np.unique(y_test)) >= 2:
            fold_qwk = cohen_kappa_score(y_test, meta_preds, weights="quadratic")
            log.info("  NAS ordinal fold %s: QWK = %.4f (n=%d)",
                     fold_ds, fold_qwk, len(y_test))

        for sid, yt, yp in zip(test_ids, y_test, meta_preds):
            all_preds.append({
                "sample_id": sid,
                "fold": fold_ds,
                "y_true_nas_group4": int(yt),
                "y_pred_nas_group4": int(yp),
                "target": "nas_group4",
            })

    if not all_preds:
        return None

    return pd.DataFrame(all_preds)


# ===================================================================
# 5. DeLong test for AUROC comparison
# ===================================================================

def delong_test(y_true, y_score1, y_score2):
    """
    DeLong test for comparing two AUROCs on the same samples.

    Implements the fast DeLong algorithm (Sun & Xu, 2014).
    Returns: z_statistic, p_value.
    """
    n1 = np.sum(y_true == 1)
    n0 = np.sum(y_true == 0)

    if n1 == 0 or n0 == 0:
        return np.nan, np.nan

    # Structural components (placement values)
    pos_mask = y_true == 1
    neg_mask = y_true == 0

    def placement_values(scores):
        """Compute placement values for DeLong test."""
        pos_scores = scores[pos_mask]
        neg_scores = scores[neg_mask]

        # V10: for each positive, fraction of negatives scored below
        v10 = np.array([
            np.mean(neg_scores < ps) + 0.5 * np.mean(neg_scores == ps)
            for ps in pos_scores
        ])
        # V01: for each negative, fraction of positives scored above
        v01 = np.array([
            np.mean(pos_scores > ns) + 0.5 * np.mean(pos_scores == ns)
            for ns in neg_scores
        ])
        return v10, v01

    v10_1, v01_1 = placement_values(y_score1)
    v10_2, v01_2 = placement_values(y_score2)

    # AUC estimates
    auc1 = np.mean(v10_1)
    auc2 = np.mean(v10_2)

    # Covariance matrix of the two AUCs
    s10_11 = np.cov(v10_1, v10_2)[0, 1] if len(v10_1) > 1 else 0
    s01_11 = np.cov(v01_1, v01_2)[0, 1] if len(v01_1) > 1 else 0
    s10_1 = np.var(v10_1, ddof=1) if len(v10_1) > 1 else 0
    s10_2 = np.var(v10_2, ddof=1) if len(v10_2) > 1 else 0
    s01_1 = np.var(v01_1, ddof=1) if len(v01_1) > 1 else 0
    s01_2 = np.var(v01_2, ddof=1) if len(v01_2) > 1 else 0

    # Variance of AUC1 - AUC2
    var_diff = (
        s10_1 / n1 + s01_1 / n0
        + s10_2 / n1 + s01_2 / n0
        - 2 * s10_11 / n1 - 2 * s01_11 / n0
    )

    if var_diff <= 0:
        return np.nan, np.nan

    z = (auc1 - auc2) / np.sqrt(var_diff)
    p = 2 * stats.norm.sf(np.abs(z))  # two-sided

    return float(z), float(p)


# ===================================================================
# 6. Comparison with single models
# ===================================================================

def build_comparison_table(stack_preds, stack_df, meta):
    """
    Compare stacking AUROC with best single model AUROC.
    Uses DeLong test for statistical significance.
    """
    rows = []

    # --- Stacking AUROC ---
    if stack_preds is not None and len(stack_preds) > 0:
        y_true = stack_preds["y_true"].values
        p_stack = stack_preds["p_stack"].values

        if len(np.unique(y_true)) >= 2:
            stack_auroc = roc_auc_score(y_true, p_stack)
            log.info("Stacking F>=3 AUROC: %.4f (n=%d)", stack_auroc, len(y_true))

            rows.append({
                "model": "stacking_ensemble",
                "target": "fib_ge3",
                "metric": "auroc",
                "value": stack_auroc,
                "n_samples": len(y_true),
            })
        else:
            stack_auroc = np.nan
    else:
        stack_auroc = np.nan

    # --- Per-base-model AUROC from stacking matrix ---
    if stack_df is not None:
        prob_cols = {"p_v3": "v3_elastic_net", "p_vae": "vae_rf",
                     "p_cbm": "concept_bottleneck", "p_mtl": "multitask_corn"}

        for col, model_name in prob_cols.items():
            if col not in stack_df.columns:
                continue
            valid = stack_df[["y_true", col]].dropna()
            if len(valid) == 0 or len(np.unique(valid["y_true"])) < 2:
                continue

            auroc = roc_auc_score(valid["y_true"].values, valid[col].values)
            rows.append({
                "model": model_name,
                "target": "fib_ge3",
                "metric": "auroc",
                "value": auroc,
                "n_samples": len(valid),
            })
            log.info("  %s F>=3 AUROC: %.4f (n=%d)", model_name, auroc, len(valid))

            # DeLong test vs stacking
            if not np.isnan(stack_auroc) and stack_preds is not None:
                # Align samples
                merged = stack_preds[["sample_id", "y_true", "p_stack"]].merge(
                    stack_df[["sample_id", col]], on="sample_id"
                ).dropna()

                if len(merged) > 10 and len(np.unique(merged["y_true"])) >= 2:
                    z_stat, p_val = delong_test(
                        merged["y_true"].values,
                        merged["p_stack"].values,
                        merged[col].values
                    )
                    rows.append({
                        "model": f"stacking_vs_{model_name}",
                        "target": "fib_ge3",
                        "metric": "delong_z",
                        "value": z_stat,
                        "n_samples": len(merged),
                    })
                    rows.append({
                        "model": f"stacking_vs_{model_name}",
                        "target": "fib_ge3",
                        "metric": "delong_p",
                        "value": p_val,
                        "n_samples": len(merged),
                    })

    # --- Load summary files for reference single-model AUROCs ---
    summary_files = {
        "v3_summary": V3_SUMMARY,
        "vae_summary": VAE_SUMMARY,
        "cbm_summary": CBM_SUMMARY,
        "mtl_summary": MTL_SUMMARY,
    }

    for name, path in summary_files.items():
        if not os.path.exists(path):
            continue
        try:
            sdf = pd.read_csv(path)
            # Extract AUROC for fib_ge3 or equivalent
            if "mean_auroc" in sdf.columns:
                for _, row in sdf.iterrows():
                    target_key = row.get("target_key", row.get("target", ""))
                    if "fib" in str(target_key).lower() or "disease" in str(target_key).lower():
                        rows.append({
                            "model": f"{name}_official",
                            "target": str(target_key),
                            "metric": "mean_auroc_official",
                            "value": row["mean_auroc"],
                            "n_samples": row.get("n_samples", np.nan),
                        })
            elif "auroc_mean" in sdf.columns:
                for _, row in sdf.iterrows():
                    task = row.get("task", "")
                    if "disease" in str(task).lower():
                        rows.append({
                            "model": f"{name}_official",
                            "target": str(task),
                            "metric": "mean_auroc_official",
                            "value": row["auroc_mean"],
                            "n_samples": np.nan,
                        })
        except Exception as e:
            log.warning("Could not parse %s: %s", name, e)

    return pd.DataFrame(rows) if rows else None


# ===================================================================
# 7. Main
# ===================================================================

def main():
    t0 = time.time()

    # Load all inputs
    log.info("Loading base model results...")
    v3_df = load_v3_predictions()
    vae_df = load_vae_classical_predictions()
    cbm_df = load_cbm_predictions()
    mtl_df = load_mtl_predictions()
    meta = load_metadata()

    # --- F>=3 Binary Stacking ---
    log.info("\n" + "=" * 60)
    log.info("PART 1: F>=3 Binary Stacking")
    log.info("=" * 60)

    stack_df = None
    stack_preds = None

    if v3_df is not None:
        stack_df = build_fib_ge3_stacking_matrix(v3_df, meta)

        if stack_df is not None:
            stack_preds = run_fib_ge3_stacking(stack_df)

    # --- NAS Ordinal Stacking ---
    log.info("\n" + "=" * 60)
    log.info("PART 2: NAS Ordinal Stacking")
    log.info("=" * 60)

    nas_preds = build_nas_ordinal_stacking_matrix(meta)

    # --- Combine results ---
    results_parts = []
    if stack_preds is not None:
        stack_preds["target"] = "fib_ge3"
        results_parts.append(stack_preds)

    if nas_preds is not None:
        results_parts.append(nas_preds)

    if results_parts:
        all_results = pd.concat(results_parts, ignore_index=True)
        all_results.to_csv(OUT_RESULTS, index=False)
        log.info("Saved stacking results: %s (%d rows)", OUT_RESULTS, len(all_results))
    else:
        log.warning("No stacking results generated")
        all_results = pd.DataFrame()

    # --- Comparison ---
    log.info("\n" + "=" * 60)
    log.info("PART 3: Stacking vs Single Model Comparison")
    log.info("=" * 60)

    comparison = build_comparison_table(stack_preds, stack_df, meta)
    if comparison is not None:
        comparison.to_csv(OUT_COMPARISON, index=False)
        log.info("Saved comparison: %s (%d rows)", OUT_COMPARISON, len(comparison))

        # Print summary
        log.info("\n--- Comparison Summary ---")
        for _, row in comparison.iterrows():
            log.info("  %-30s  %s = %.4f",
                     row["model"], row["metric"], row["value"])
    else:
        log.warning("No comparison table generated")

    # --- NAS Ordinal QWK Summary ---
    if nas_preds is not None and len(nas_preds) > 0:
        y_true_ord = nas_preds["y_true_nas_group4"].values
        y_pred_ord = nas_preds["y_pred_nas_group4"].values
        overall_qwk = cohen_kappa_score(y_true_ord, y_pred_ord, weights="quadratic")
        log.info("\nNAS 4-group ordinal stacking QWK: %.4f (n=%d)",
                 overall_qwk, len(nas_preds))

    elapsed = time.time() - t0
    log.info("\n" + "=" * 70)
    log.info("108_stacking_ensemble.py completed in %.1f seconds", elapsed)
    log.info("=" * 70)


if __name__ == "__main__":
    main()
