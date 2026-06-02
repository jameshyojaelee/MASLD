#!/usr/bin/env python3
"""
108c_proper_stacking.py
Proper stacking ensemble with DIVERSE feature representations.

CRITICAL CONTEXT:
  108  had data leakage (full-model embeddings used in LOCO-CV).
  108b fixed leakage but stacked 3 classifiers on the SAME 64-dim embeddings
       → no diversity → QWK 0.343 (worse than single models).

KEY INSIGHT: Stacking gains come from combining base models with DIFFERENT
inductive biases operating on DIFFERENT feature representations.

Base models (per LOCO fold):
  1. Expression elastic net  — top 3K genes by training-fold variance, rank-transformed
  2. NAS-VAE RF              — LOCO-specific VAE embeddings (64-dim)
  3. Pathway elastic net     — 1791 ssGSEA pathway scores (unsupervised → no leakage)
  4. TF activity XGBoost     — 295 decoupleR TF activity scores (unsupervised → no leakage)

Each base model produces out-of-fold predicted class probabilities.
Meta-learner: logistic regression on stacked base-model probabilities,
trained on OTHER folds' stacked predictions via nested LOCO.

Targets:
  - NAS 4-group ordinal (QWK)
  - NAS>=5 binary (AUROC)
  - F>=3 binary (AUROC)

ANTI-LEAKAGE RULES:
  - LOCO-specific NAS embeddings from nas_embeddings_loco/ only
  - Expression features selected within each fold (training-fold variance)
  - Pathway scores (ssGSEA) and TF activities (decoupleR) are unsupervised — OK to use full matrix
  - Meta-learner trained ONLY on properly held-out base model predictions

Output:
  - proper_stacking_results.csv      (per-sample predictions for all targets)
  - proper_stacking_summary.csv      (aggregated metrics)

SLURM: cpu, 8 CPUs, 64G, 48h

Usage:
  sbatch --job-name=stg108c_stacking \
         --partition=cpu --cpus-per-task=8 --mem=64G --time=48:00:00 \
         --output=logs/108c_proper_stacking_%j.out \
         --error=logs/108c_proper_stacking_%j.err \
         --wrap="bash -c 'eval \"\\$(micromamba shell hook --shell bash)\" && \\
                 micromamba activate rapids_singlecell && \\
                 cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts && \\
                 python 108c_proper_stacking.py'"
"""

# ============================================================
# Environment guards — prevent CUDA crash on CPU-only nodes
# ============================================================
import os
os.environ["RAPIDS_NO_INITIALIZE"] = "1"
os.environ["CUDF_PANDAS"] = "0"

import sys
import time
import warnings
import subprocess
import tempfile
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, confusion_matrix

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ============================================================
# Reproducibility
# ============================================================
SEED = 42
np.random.seed(SEED)

# ============================================================
# Paths
# ============================================================
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

RSCRIPT = "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"

# Input files
META_PATH = os.path.join(OUTDIR, "modeling_metadata.csv")
NAS_EMB_DIR = os.path.join(OUTDIR, "nas_embeddings_loco")
FIB_EMB_DIR = os.path.join(OUTDIR, "embeddings_loco")
SSGSEA_CSV = os.path.join(OUTDIR, "pathway_scores_ssgsea.csv")
SSGSEA_RDS = os.path.join(OUTDIR, "pathway_scores_ssgsea.rds")
TF_CSV = os.path.join(OUTDIR, "tf_activity_features.csv")
TF_RDS = os.path.join(OUTDIR, "tf_activity_features.rds")
PATHWAY_NAMES = os.path.join(OUTDIR, "pathway_feature_names.csv")
TF_NAMES = os.path.join(OUTDIR, "tf_activity_names.csv")
RANK_EXPR_FULL_RDS = os.path.join(OUTDIR, "rank_expression_full.rds")

# Outputs
OUT_RESULTS = os.path.join(OUTDIR, "proper_stacking_results.csv")
OUT_SUMMARY = os.path.join(OUTDIR, "proper_stacking_summary.csv")

# ============================================================
# Logging
# ============================================================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(os.path.join(LOGDIR, "108c_proper_stacking.log")),
    ],
)
log = logging.getLogger(__name__)

t0 = time.time()
log.info("=== 108c: Proper Stacking Ensemble (Diverse Features) ===")


# ============================================================
# Helper: QWK
# ============================================================
def quadratic_weighted_kappa(y_true, y_pred, n_classes):
    """Compute quadratic weighted kappa (ordinal agreement)."""
    cm = confusion_matrix(y_true, y_pred, labels=list(range(n_classes)))
    n = cm.sum()
    if n == 0:
        return 0.0
    w = np.zeros((n_classes, n_classes))
    for i in range(n_classes):
        for j in range(n_classes):
            w[i, j] = (i - j) ** 2 / (n_classes - 1) ** 2
    expected = np.outer(cm.sum(axis=1), cm.sum(axis=0)) / n
    num = (w * cm).sum()
    den = (w * expected).sum()
    return 1.0 - num / den if den > 0 else 1.0


# ============================================================
# Helper: Load RDS via Rscript subprocess → CSV
# ============================================================
def load_rds_matrix(rds_path, transpose=False):
    """Load an RDS matrix file and return as pandas DataFrame.

    Parameters
    ----------
    rds_path : str
        Path to the .rds file containing a matrix/data.frame.
    transpose : bool
        If True, transpose the matrix before writing CSV
        (e.g., for ssGSEA: pathways x samples → samples x pathways).

    Returns
    -------
    pd.DataFrame with row names as index, column names as columns.
    """
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
        tmp_csv = tmp.name

    t_flag = "TRUE" if transpose else "FALSE"
    r_code = f"""
x <- readRDS("{rds_path}")
if ({t_flag}) x <- t(x)
write.csv(x, "{tmp_csv}", row.names = TRUE)
"""
    result = subprocess.run(
        [RSCRIPT, "-e", r_code],
        capture_output=True, text=True, timeout=300,
    )
    if result.returncode != 0:
        log.error(f"Rscript failed for {rds_path}: {result.stderr}")
        os.unlink(tmp_csv)
        raise RuntimeError(f"Failed to load RDS: {rds_path}")

    df = pd.read_csv(tmp_csv, index_col=0)
    os.unlink(tmp_csv)
    return df


# ============================================================
# Load metadata
# ============================================================
log.info("Loading metadata...")
meta = pd.read_csv(META_PATH)
meta["sample_id"] = meta["sample_id"].astype(str).str.strip()

# NAS LOCO folds (5 cohorts with NAS scores)
NAS_DATASETS = ["GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066"]
nas_meta = meta[
    (meta["loco_fold_nas"] != "excluded") & (meta["nas_group4"] >= 0)
].copy()
log.info(f"NAS-eligible samples: {len(nas_meta)}, folds: {nas_meta['loco_fold_nas'].nunique()}")

# Fibrosis LOCO folds (6 cohorts with fibrosis staging)
FIB_DATASETS = ["GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066", "GSE240729"]
fib_meta = meta[
    (meta["loco_fold_fibrosis"] != "excluded") & (meta["fib_ge3"] >= 0)
].copy()
log.info(f"Fibrosis-eligible samples: {len(fib_meta)}, folds: {fib_meta['loco_fold_fibrosis'].nunique()}")

# ============================================================
# Load feature matrices
# ============================================================

# --- 1. Full logCPM expression (for within-fold feature selection) ---
log.info("Loading full expression matrix (34,453 genes x 1,444 samples) from RDS...")
expr_full = load_rds_matrix(RANK_EXPR_FULL_RDS, transpose=False)
# This is genes x samples; transpose to samples x genes for sklearn
expr_full = expr_full.T
expr_full.index = expr_full.index.astype(str).str.strip()
log.info(f"  Expression: {expr_full.shape[0]} samples x {expr_full.shape[1]} genes")

# --- 2. NAS-VAE LOCO embeddings ---
log.info("Loading NAS-VAE LOCO embeddings...")
nas_loco_embeddings = {}
for fold in NAS_DATASETS:
    emb_file = os.path.join(NAS_EMB_DIR, f"nas_embeddings_loco_{fold}.csv")
    if os.path.exists(emb_file):
        df = pd.read_csv(emb_file)
        df["sample_id"] = df["sample_id"].astype(str).str.strip()
        df = df.set_index("sample_id")
        nas_loco_embeddings[fold] = df
        log.info(f"  NAS LOCO {fold}: {df.shape}")
    else:
        log.warning(f"  Missing NAS LOCO embeddings for {fold}")

# --- 2b. General LOCO embeddings (for fibrosis folds) ---
log.info("Loading general LOCO embeddings for fibrosis folds...")
fib_loco_embeddings = {}
for fold in FIB_DATASETS:
    emb_file = os.path.join(FIB_EMB_DIR, f"embeddings_loco_{fold}.csv")
    if os.path.exists(emb_file):
        df = pd.read_csv(emb_file)
        df["sample_id"] = df["sample_id"].astype(str).str.strip()
        df = df.set_index("sample_id")
        fib_loco_embeddings[fold] = df
        log.info(f"  Fibrosis LOCO {fold}: {df.shape}")
    else:
        log.warning(f"  Missing fibrosis LOCO embeddings for {fold}")

# --- 3. ssGSEA pathway scores ---
log.info("Loading ssGSEA pathway scores...")
# The CSV is pathways (rows) x samples (columns) with no row-name column.
# Use the RDS which preserves row/col names properly, and transpose.
if os.path.exists(SSGSEA_RDS):
    pathway_df = load_rds_matrix(SSGSEA_RDS, transpose=True)
    pathway_df.index = pathway_df.index.astype(str).str.strip()
    log.info(f"  ssGSEA (from RDS, transposed): {pathway_df.shape[0]} samples x {pathway_df.shape[1]} pathways")
else:
    # Fallback: CSV with sample IDs as columns, pathway names from external file
    log.info("  Falling back to CSV for ssGSEA...")
    raw = pd.read_csv(SSGSEA_CSV, header=0)
    # Columns are sample IDs; rows are pathways. Transpose.
    raw.columns = [str(c).strip().strip('"') for c in raw.columns]
    pnames = pd.read_csv(PATHWAY_NAMES)
    raw.index = pnames["pathway"].values[:len(raw)]
    pathway_df = raw.T
    pathway_df.index = pathway_df.index.astype(str).str.strip()
    log.info(f"  ssGSEA (from CSV, transposed): {pathway_df.shape[0]} samples x {pathway_df.shape[1]} pathways")

# --- 4. TF activity scores ---
log.info("Loading TF activity scores...")
if os.path.exists(TF_RDS):
    tf_df = load_rds_matrix(TF_RDS, transpose=False)
    tf_df.index = tf_df.index.astype(str).str.strip()
    log.info(f"  TF activity (from RDS): {tf_df.shape[0]} samples x {tf_df.shape[1]} TFs")
else:
    # Fallback: CSV with TF names as columns, samples as rows, no sample_id column
    log.info("  Falling back to CSV for TF activity...")
    raw = pd.read_csv(TF_CSV, header=0)
    raw.columns = [str(c).strip().strip('"') for c in raw.columns]
    # Align to metadata sample order
    tf_df = raw.copy()
    # Use metadata sample_id as index (same order)
    tf_df.index = meta["sample_id"].values[:len(tf_df)]
    tf_df.index = tf_df.index.astype(str).str.strip()
    log.info(f"  TF activity (from CSV): {tf_df.shape[0]} samples x {tf_df.shape[1]} TFs")

# ============================================================
# Verify XGBoost availability
# ============================================================
try:
    from xgboost import XGBClassifier
    HAS_XGB = True
    log.info("XGBoost available")
except ImportError:
    from sklearn.ensemble import GradientBoostingClassifier
    HAS_XGB = False
    log.warning("XGBoost not available — using sklearn GradientBoostingClassifier as fallback")

# ============================================================
# Helper: train base model and get predicted probabilities
# ============================================================
K_FEATURES = 3000


def train_expression_elastic_net(train_ids, test_ids, y_train, n_classes, is_binary=False):
    """Base model 1: Expression elastic net with within-fold feature selection.

    1. Select top 3K genes by variance on training samples only.
    2. Rank-transform expression values.
    3. Fit elastic net (logistic for binary, multinomial for multi-class).
    4. Return predicted probabilities on test set.
    """
    # Filter to samples present in expression matrix
    avail_train = [s for s in train_ids if s in expr_full.index]
    avail_test = [s for s in test_ids if s in expr_full.index]

    if len(avail_test) < 2 or len(avail_train) < 10:
        return None, avail_test

    X_train_raw = expr_full.loc[avail_train]
    X_test_raw = expr_full.loc[avail_test]

    # Step 1: top 3K genes by TRAINING variance
    gene_var = X_train_raw.var(axis=0)
    top_genes = gene_var.nlargest(K_FEATURES).index.tolist()

    X_train_sel = X_train_raw[top_genes]
    X_test_sel = X_test_raw[top_genes]

    # Step 2: rank transform per sample (within each sample, rank genes)
    X_train_rank = X_train_sel.rank(axis=1, method="average")
    X_train_rank = X_train_rank / X_train_rank.shape[1]
    X_test_rank = X_test_sel.rank(axis=1, method="average")
    X_test_rank = X_test_rank / X_test_rank.shape[1]

    # Step 3: scale
    scaler = StandardScaler()
    X_tr = scaler.fit_transform(X_train_rank.values)
    X_te = scaler.transform(X_test_rank.values)

    # Align labels to available samples
    label_map = dict(zip(train_ids, y_train))
    y_tr = np.array([label_map[s] for s in avail_train])

    # Step 4: elastic net
    if is_binary:
        model = LogisticRegression(
            penalty="elasticnet", solver="saga", l1_ratio=0.5,
            C=1.0, max_iter=10000, class_weight="balanced",
            random_state=SEED,
        )
    else:
        model = LogisticRegression(
            penalty="elasticnet", solver="saga", l1_ratio=0.5,
            C=1.0, max_iter=10000,
            class_weight="balanced", random_state=SEED,
        )

    model.fit(X_tr, y_tr)
    probs = model.predict_proba(X_te)
    return probs, avail_test


def train_vae_rf(train_ids, test_ids, y_train, fold_emb, n_classes):
    """Base model 2: Random forest on LOCO-specific VAE embeddings."""
    emb_cols = [c for c in fold_emb.columns if c.startswith("z")]

    avail_train = [s for s in train_ids if s in fold_emb.index]
    avail_test = [s for s in test_ids if s in fold_emb.index]

    if len(avail_test) < 2 or len(avail_train) < 10:
        return None, avail_test

    X_tr = fold_emb.loc[avail_train, emb_cols].values.astype(np.float32)
    X_te = fold_emb.loc[avail_test, emb_cols].values.astype(np.float32)

    label_map = dict(zip(train_ids, y_train))
    y_tr = np.array([label_map[s] for s in avail_train])

    scaler = StandardScaler()
    X_tr_s = scaler.fit_transform(X_tr)
    X_te_s = scaler.transform(X_te)

    model = RandomForestClassifier(
        n_estimators=500, class_weight="balanced",
        random_state=SEED, n_jobs=-1, max_features="sqrt",
    )
    model.fit(X_tr_s, y_tr)
    probs = model.predict_proba(X_te_s)
    return probs, avail_test


def train_pathway_elastic_net(train_ids, test_ids, y_train, n_classes, is_binary=False):
    """Base model 3: Elastic net on ssGSEA pathway scores."""
    avail_train = [s for s in train_ids if s in pathway_df.index]
    avail_test = [s for s in test_ids if s in pathway_df.index]

    if len(avail_test) < 2 or len(avail_train) < 10:
        return None, avail_test

    X_tr = pathway_df.loc[avail_train].values.astype(np.float32)
    X_te = pathway_df.loc[avail_test].values.astype(np.float32)

    label_map = dict(zip(train_ids, y_train))
    y_tr = np.array([label_map[s] for s in avail_train])

    scaler = StandardScaler()
    X_tr_s = scaler.fit_transform(X_tr)
    X_te_s = scaler.transform(X_te)

    if is_binary:
        model = LogisticRegression(
            penalty="elasticnet", solver="saga", l1_ratio=0.5,
            C=1.0, max_iter=10000, class_weight="balanced",
            random_state=SEED,
        )
    else:
        model = LogisticRegression(
            penalty="elasticnet", solver="saga", l1_ratio=0.5,
            C=1.0, max_iter=10000,
            class_weight="balanced", random_state=SEED,
        )

    model.fit(X_tr_s, y_tr)
    probs = model.predict_proba(X_te_s)
    return probs, avail_test


def train_tf_xgboost(train_ids, test_ids, y_train, n_classes, is_binary=False):
    """Base model 4: XGBoost (or GradientBoosting fallback) on TF activity scores."""
    avail_train = [s for s in train_ids if s in tf_df.index]
    avail_test = [s for s in test_ids if s in tf_df.index]

    if len(avail_test) < 2 or len(avail_train) < 10:
        return None, avail_test

    X_tr = tf_df.loc[avail_train].values.astype(np.float32)
    X_te = tf_df.loc[avail_test].values.astype(np.float32)

    label_map = dict(zip(train_ids, y_train))
    y_tr = np.array([label_map[s] for s in avail_train])

    scaler = StandardScaler()
    X_tr_s = scaler.fit_transform(X_tr)
    X_te_s = scaler.transform(X_te)

    if HAS_XGB:
        if is_binary:
            model = XGBClassifier(
                n_estimators=300, max_depth=4, learning_rate=0.05,
                subsample=0.8, colsample_bytree=0.8,
                objective="binary:logistic", eval_metric="logloss",
                random_state=SEED, n_jobs=-1, verbosity=0,
            )
        else:
            model = XGBClassifier(
                n_estimators=300, max_depth=4, learning_rate=0.05,
                subsample=0.8, colsample_bytree=0.8,
                objective="multi:softprob", eval_metric="mlogloss",
                num_class=n_classes,
                random_state=SEED, n_jobs=-1, verbosity=0,
            )
    else:
        model = GradientBoostingClassifier(
            n_estimators=300, max_depth=4, learning_rate=0.05,
            subsample=0.8, random_state=SEED,
        )

    model.fit(X_tr_s, y_tr)
    probs = model.predict_proba(X_te_s)
    return probs, avail_test


# ============================================================
# Run stacking for one target (generic function)
# ============================================================
def run_stacking(target_name, eligible_meta, fold_col, label_col, folds,
                 loco_embeddings, n_classes, is_binary=False):
    """Run full LOCO stacking pipeline for one target.

    For each outer LOCO fold:
      1. Train 4 base models on training folds → predict held-out fold probs
      2. Collect all out-of-fold base model predictions

    Then nested meta-learner:
      For each outer fold, train meta-learner on the OTHER folds' stacked
      predictions and predict this fold → final out-of-fold meta-predictions.

    Returns per-sample DataFrame with predictions.
    """
    log.info(f"\n{'='*60}")
    log.info(f"Target: {target_name} ({'binary' if is_binary else f'{n_classes}-class ordinal'})")
    log.info(f"{'='*60}")

    label_map = dict(zip(eligible_meta["sample_id"], eligible_meta[label_col]))
    fold_map = dict(zip(eligible_meta["sample_id"], eligible_meta[fold_col]))

    base_model_names = ["expr_enet", "vae_rf", "pathway_enet", "tf_xgb"]
    n_base = len(base_model_names)

    # Phase 1: Collect out-of-fold base model predictions for ALL folds
    # For each fold, each base model predicts the held-out fold
    # Store: {sample_id: {model_name: probability_vector}}
    oof_base_preds = {}   # sample_id → np.array of shape (n_base * n_classes,)
    oof_true_labels = {}  # sample_id → int label
    oof_folds = {}        # sample_id → fold name
    per_fold_base_metrics = []

    for fold in folds:
        if fold not in loco_embeddings:
            log.warning(f"  Skipping {fold}: no LOCO embeddings")
            continue

        # Partition into train/test
        test_ids = eligible_meta[eligible_meta[fold_col] == fold]["sample_id"].values
        train_ids = eligible_meta[eligible_meta[fold_col] != fold]["sample_id"].values
        y_train = np.array([label_map[s] for s in train_ids])

        if len(test_ids) < 5 or len(train_ids) < 20:
            log.warning(f"  Skipping {fold}: too few samples (test={len(test_ids)}, train={len(train_ids)})")
            continue

        log.info(f"  Fold {fold}: train={len(train_ids)}, test={len(test_ids)}")

        fold_emb = loco_embeddings[fold]

        # Train each base model and collect out-of-fold predictions
        base_results = {}

        # Base model 1: Expression elastic net
        probs, avail = train_expression_elastic_net(
            train_ids, test_ids, y_train, n_classes, is_binary=is_binary
        )
        if probs is not None:
            base_results["expr_enet"] = (probs, avail)
            pred = probs.argmax(axis=1) if not is_binary else (probs[:, 1] > 0.5).astype(int)
            if not is_binary:
                qwk = quadratic_weighted_kappa(
                    np.array([label_map[s] for s in avail]), pred, n_classes
                )
                log.info(f"    expr_enet: QWK={qwk:.3f}")
            else:
                y_te = np.array([label_map[s] for s in avail])
                if len(np.unique(y_te)) == 2:
                    auc = roc_auc_score(y_te, probs[:, 1])
                    log.info(f"    expr_enet: AUROC={auc:.3f}")
        else:
            log.warning(f"    expr_enet: SKIPPED (insufficient samples)")

        # Base model 2: VAE RF
        probs, avail = train_vae_rf(
            train_ids, test_ids, y_train, fold_emb, n_classes
        )
        if probs is not None:
            base_results["vae_rf"] = (probs, avail)
            pred = probs.argmax(axis=1) if not is_binary else (probs[:, 1] > 0.5).astype(int)
            if not is_binary:
                qwk = quadratic_weighted_kappa(
                    np.array([label_map[s] for s in avail]), pred, n_classes
                )
                log.info(f"    vae_rf: QWK={qwk:.3f}")
            else:
                y_te = np.array([label_map[s] for s in avail])
                if len(np.unique(y_te)) == 2:
                    auc = roc_auc_score(y_te, probs[:, 1])
                    log.info(f"    vae_rf: AUROC={auc:.3f}")
        else:
            log.warning(f"    vae_rf: SKIPPED (insufficient samples)")

        # Base model 3: Pathway elastic net
        probs, avail = train_pathway_elastic_net(
            train_ids, test_ids, y_train, n_classes, is_binary=is_binary
        )
        if probs is not None:
            base_results["pathway_enet"] = (probs, avail)
            pred = probs.argmax(axis=1) if not is_binary else (probs[:, 1] > 0.5).astype(int)
            if not is_binary:
                qwk = quadratic_weighted_kappa(
                    np.array([label_map[s] for s in avail]), pred, n_classes
                )
                log.info(f"    pathway_enet: QWK={qwk:.3f}")
            else:
                y_te = np.array([label_map[s] for s in avail])
                if len(np.unique(y_te)) == 2:
                    auc = roc_auc_score(y_te, probs[:, 1])
                    log.info(f"    pathway_enet: AUROC={auc:.3f}")
        else:
            log.warning(f"    pathway_enet: SKIPPED (insufficient samples)")

        # Base model 4: TF XGBoost
        probs, avail = train_tf_xgboost(
            train_ids, test_ids, y_train, n_classes, is_binary=is_binary
        )
        if probs is not None:
            base_results["tf_xgb"] = (probs, avail)
            pred = probs.argmax(axis=1) if not is_binary else (probs[:, 1] > 0.5).astype(int)
            if not is_binary:
                qwk = quadratic_weighted_kappa(
                    np.array([label_map[s] for s in avail]), pred, n_classes
                )
                log.info(f"    tf_xgb: QWK={qwk:.3f}")
            else:
                y_te = np.array([label_map[s] for s in avail])
                if len(np.unique(y_te)) == 2:
                    auc = roc_auc_score(y_te, probs[:, 1])
                    log.info(f"    tf_xgb: AUROC={auc:.3f}")
        else:
            log.warning(f"    tf_xgb: SKIPPED (insufficient samples)")

        # Find samples present in ALL base models for stacking
        if len(base_results) < 2:
            log.warning(f"  Skipping {fold}: fewer than 2 base models produced predictions")
            continue

        # Find common test samples across all available base models
        available_models = list(base_results.keys())
        sample_sets = [set(base_results[m][1]) for m in available_models]
        common_samples = list(sorted(set.intersection(*sample_sets)))

        if len(common_samples) < 2:
            log.warning(f"  Skipping {fold}: only {len(common_samples)} common test samples")
            continue

        log.info(f"    Stacking: {len(common_samples)} common samples, {len(available_models)} base models")

        # Build stacked feature matrix for these common samples
        for sid in common_samples:
            prob_vec = []
            for mname in base_model_names:
                if mname in base_results:
                    probs_arr, avail_arr = base_results[mname]
                    idx = list(avail_arr).index(sid)
                    prob_vec.extend(probs_arr[idx].tolist())
                else:
                    # Pad with uniform probabilities for missing models
                    prob_vec.extend([1.0 / n_classes] * n_classes)

            oof_base_preds[sid] = np.array(prob_vec)
            oof_true_labels[sid] = label_map[sid]
            oof_folds[sid] = fold

    # Phase 2: Nested meta-learner
    # For each fold, train meta-learner on stacked predictions from OTHER folds
    log.info(f"\n  Phase 2: Training meta-learner (nested LOCO)...")

    all_stacked_samples = list(oof_base_preds.keys())
    n_stacked_features = len(next(iter(oof_base_preds.values()))) if all_stacked_samples else 0

    if len(all_stacked_samples) < 20:
        log.error(f"  Only {len(all_stacked_samples)} samples with stacked predictions — aborting target {target_name}")
        return pd.DataFrame()

    # Build full stacked matrix
    X_stacked = np.array([oof_base_preds[s] for s in all_stacked_samples])
    y_stacked = np.array([oof_true_labels[s] for s in all_stacked_samples])
    fold_stacked = np.array([oof_folds[s] for s in all_stacked_samples])

    log.info(f"  Stacked matrix: {X_stacked.shape[0]} samples x {X_stacked.shape[1]} features")
    log.info(f"  Label distribution: {dict(zip(*np.unique(y_stacked, return_counts=True)))}")

    # Meta-learner: nested LOCO (train on all-but-one fold, predict held-out)
    meta_preds = np.full(len(all_stacked_samples), np.nan)
    meta_probs = np.full(len(all_stacked_samples), np.nan)
    meta_probs_full = np.full((len(all_stacked_samples), n_classes), np.nan)

    unique_folds = sorted(set(fold_stacked))

    for meta_fold in unique_folds:
        meta_test_mask = fold_stacked == meta_fold
        meta_train_mask = ~meta_test_mask

        X_meta_train = X_stacked[meta_train_mask]
        y_meta_train = y_stacked[meta_train_mask]
        X_meta_test = X_stacked[meta_test_mask]

        if X_meta_test.shape[0] < 2 or len(np.unique(y_meta_train)) < 2:
            log.warning(f"    Meta-fold {meta_fold}: skipped (insufficient samples or classes)")
            continue

        # Scale stacked features
        meta_scaler = StandardScaler()
        X_mt_s = meta_scaler.fit_transform(X_meta_train)
        X_mte_s = meta_scaler.transform(X_meta_test)

        if is_binary:
            meta_model = LogisticRegression(
                C=1.0, max_iter=5000, class_weight="balanced",
                random_state=SEED,
            )
        else:
            meta_model = LogisticRegression(
                C=1.0, max_iter=5000,
                class_weight="balanced", random_state=SEED,
            )

        meta_model.fit(X_mt_s, y_meta_train)

        probs_meta = meta_model.predict_proba(X_mte_s)
        preds_meta = meta_model.predict(X_mte_s)

        meta_preds[meta_test_mask] = preds_meta
        meta_probs_full[meta_test_mask, :probs_meta.shape[1]] = probs_meta
        if is_binary and probs_meta.shape[1] == 2:
            meta_probs[meta_test_mask] = probs_meta[:, 1]

        if not is_binary:
            fold_qwk = quadratic_weighted_kappa(
                y_stacked[meta_test_mask], preds_meta, n_classes
            )
            log.info(f"    Meta-fold {meta_fold}: QWK={fold_qwk:.3f} (n={meta_test_mask.sum()})")
        else:
            y_te = y_stacked[meta_test_mask]
            if len(np.unique(y_te)) == 2:
                fold_auc = roc_auc_score(y_te, probs_meta[:, 1])
                log.info(f"    Meta-fold {meta_fold}: AUROC={fold_auc:.3f} (n={meta_test_mask.sum()})")

    # Also compute simple average ensemble (no meta-learner) for comparison
    # Average the base model probabilities per class
    avg_probs_full = np.zeros((len(all_stacked_samples), n_classes))
    n_base_actual = n_stacked_features // n_classes
    for b in range(n_base_actual):
        start = b * n_classes
        end = start + n_classes
        avg_probs_full += X_stacked[:, start:end]
    avg_probs_full /= n_base_actual
    avg_preds = avg_probs_full.argmax(axis=1)

    # Phase 3: Compute overall metrics
    valid_mask = ~np.isnan(meta_preds)
    results_rows = []

    if valid_mask.sum() > 0:
        y_valid = y_stacked[valid_mask]
        pred_valid = meta_preds[valid_mask].astype(int)

        if not is_binary:
            overall_qwk_meta = quadratic_weighted_kappa(y_valid, pred_valid, n_classes)
            overall_qwk_avg = quadratic_weighted_kappa(y_valid, avg_preds[valid_mask], n_classes)
            log.info(f"\n  *** {target_name} META-LEARNER QWK: {overall_qwk_meta:.4f} ***")
            log.info(f"  *** {target_name} AVERAGE ENSEMBLE QWK: {overall_qwk_avg:.4f} ***")
            log.info(f"      ({valid_mask.sum()} samples, {len(unique_folds)} folds)")
        else:
            if len(np.unique(y_valid)) == 2:
                probs_valid = meta_probs[valid_mask]
                overall_auc_meta = roc_auc_score(y_valid, probs_valid)
                avg_prob_binary = avg_probs_full[valid_mask, 1] if n_classes == 2 else avg_probs_full[valid_mask, -1]
                overall_auc_avg = roc_auc_score(y_valid, avg_prob_binary)
                log.info(f"\n  *** {target_name} META-LEARNER AUROC: {overall_auc_meta:.4f} ***")
                log.info(f"  *** {target_name} AVERAGE ENSEMBLE AUROC: {overall_auc_avg:.4f} ***")
                log.info(f"      ({valid_mask.sum()} samples, {len(unique_folds)} folds)")

        # Build per-sample results
        for i, sid in enumerate(all_stacked_samples):
            if not valid_mask[i]:
                continue
            row = {
                "sample_id": sid,
                "target": target_name,
                "fold": oof_folds[sid],
                "true_label": int(oof_true_labels[sid]),
                "meta_pred": int(meta_preds[i]),
                "avg_pred": int(avg_preds[i]),
            }
            # Add meta-learner probabilities
            for c in range(n_classes):
                if not np.isnan(meta_probs_full[i, c]):
                    row[f"meta_prob_c{c}"] = float(meta_probs_full[i, c])
            # Add average ensemble probabilities
            for c in range(n_classes):
                row[f"avg_prob_c{c}"] = float(avg_probs_full[i, c])
            results_rows.append(row)

    return pd.DataFrame(results_rows)


# ============================================================
# Target 1: NAS 4-group ordinal
# ============================================================
nas4_results = run_stacking(
    target_name="nas_4group",
    eligible_meta=nas_meta,
    fold_col="loco_fold_nas",
    label_col="nas_group4",
    folds=NAS_DATASETS,
    loco_embeddings=nas_loco_embeddings,
    n_classes=4,
    is_binary=False,
)

# ============================================================
# Target 2: NAS>=5 binary
# ============================================================
nas_ge5_meta = nas_meta[nas_meta["nas_ge5"] >= 0].copy()
nas_ge5_results = run_stacking(
    target_name="nas_ge5",
    eligible_meta=nas_ge5_meta,
    fold_col="loco_fold_nas",
    label_col="nas_ge5",
    folds=NAS_DATASETS,
    loco_embeddings=nas_loco_embeddings,
    n_classes=2,
    is_binary=True,
)

# ============================================================
# Target 3: F>=3 binary
# ============================================================
fib_ge3_meta = fib_meta[fib_meta["fib_ge3"] >= 0].copy()
fib_ge3_results = run_stacking(
    target_name="fib_ge3",
    eligible_meta=fib_ge3_meta,
    fold_col="loco_fold_fibrosis",
    label_col="fib_ge3",
    folds=FIB_DATASETS,
    loco_embeddings=fib_loco_embeddings,
    n_classes=2,
    is_binary=True,
)

# ============================================================
# Combine and save results
# ============================================================
log.info(f"\n{'='*60}")
log.info("Combining results and saving...")

all_results = pd.concat(
    [df for df in [nas4_results, nas_ge5_results, fib_ge3_results] if len(df) > 0],
    ignore_index=True,
)

if len(all_results) > 0:
    all_results.to_csv(OUT_RESULTS, index=False)
    log.info(f"Saved per-sample results: {OUT_RESULTS} ({len(all_results)} rows)")
else:
    log.error("No results to save!")

# ============================================================
# Summary table
# ============================================================
summary_rows = []

for target_name, df, n_classes, is_binary in [
    ("nas_4group", nas4_results, 4, False),
    ("nas_ge5", nas_ge5_results, 2, True),
    ("fib_ge3", fib_ge3_results, 2, True),
]:
    if len(df) == 0:
        summary_rows.append({
            "target": target_name, "metric": "QWK" if not is_binary else "AUROC",
            "meta_learner": np.nan, "avg_ensemble": np.nan,
            "n_samples": 0, "n_folds": 0,
        })
        continue

    y_true = df["true_label"].values
    n_samples = len(df)
    n_folds = df["fold"].nunique()

    if not is_binary:
        # QWK
        meta_qwk = quadratic_weighted_kappa(y_true, df["meta_pred"].values, n_classes)
        avg_qwk = quadratic_weighted_kappa(y_true, df["avg_pred"].values, n_classes)
        summary_rows.append({
            "target": target_name, "metric": "QWK",
            "meta_learner": round(meta_qwk, 4),
            "avg_ensemble": round(avg_qwk, 4),
            "n_samples": n_samples, "n_folds": n_folds,
        })
    else:
        # AUROC
        prob_col = "meta_prob_c1"
        avg_prob_col = "avg_prob_c1"
        if prob_col in df.columns and avg_prob_col in df.columns:
            if len(np.unique(y_true)) == 2:
                meta_auc = roc_auc_score(y_true, df[prob_col].values)
                avg_auc = roc_auc_score(y_true, df[avg_prob_col].values)
                summary_rows.append({
                    "target": target_name, "metric": "AUROC",
                    "meta_learner": round(meta_auc, 4),
                    "avg_ensemble": round(avg_auc, 4),
                    "n_samples": n_samples, "n_folds": n_folds,
                })
            else:
                summary_rows.append({
                    "target": target_name, "metric": "AUROC",
                    "meta_learner": np.nan, "avg_ensemble": np.nan,
                    "n_samples": n_samples, "n_folds": n_folds,
                })
        else:
            summary_rows.append({
                "target": target_name, "metric": "AUROC",
                "meta_learner": np.nan, "avg_ensemble": np.nan,
                "n_samples": n_samples, "n_folds": n_folds,
            })

    # Per-fold breakdown
    for fold in sorted(df["fold"].unique()):
        fdf = df[df["fold"] == fold]
        y_f = fdf["true_label"].values
        n_f = len(fdf)

        if not is_binary:
            fold_meta_qwk = quadratic_weighted_kappa(y_f, fdf["meta_pred"].values, n_classes)
            fold_avg_qwk = quadratic_weighted_kappa(y_f, fdf["avg_pred"].values, n_classes)
            summary_rows.append({
                "target": f"{target_name}_fold_{fold}", "metric": "QWK",
                "meta_learner": round(fold_meta_qwk, 4),
                "avg_ensemble": round(fold_avg_qwk, 4),
                "n_samples": n_f, "n_folds": 1,
            })
        else:
            prob_col = "meta_prob_c1"
            avg_prob_col = "avg_prob_c1"
            if prob_col in fdf.columns and len(np.unique(y_f)) == 2:
                fold_meta_auc = roc_auc_score(y_f, fdf[prob_col].values)
                fold_avg_auc = roc_auc_score(y_f, fdf[avg_prob_col].values)
                summary_rows.append({
                    "target": f"{target_name}_fold_{fold}", "metric": "AUROC",
                    "meta_learner": round(fold_meta_auc, 4),
                    "avg_ensemble": round(fold_avg_auc, 4),
                    "n_samples": n_f, "n_folds": 1,
                })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(OUT_SUMMARY, index=False)
log.info(f"Saved summary: {OUT_SUMMARY}")

# Print final summary
log.info(f"\n{'='*60}")
log.info("FINAL SUMMARY")
log.info(f"{'='*60}")
for _, row in summary_df.iterrows():
    if "_fold_" not in str(row["target"]):
        log.info(f"  {row['target']:15s} {row['metric']:6s} | "
                 f"meta={row['meta_learner']:.4f} avg={row['avg_ensemble']:.4f} | "
                 f"n={row['n_samples']} folds={row['n_folds']}")

elapsed = time.time() - t0
log.info(f"\nTotal runtime: {elapsed/60:.1f} minutes")
log.info("=== 108c_proper_stacking.py completed ===")
