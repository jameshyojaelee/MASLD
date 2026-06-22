#!/usr/bin/env python3
"""190_multiprogram_predictor.py — Multi-task ordinal prediction via concept
bottleneck and elastic net, with modality ablation and random gene baselines.

Replaces F>=3 binary prediction with multi-target ordinal prediction across:
  - Fibrosis (F0-F4, 5 levels, N=1128, 6 LOCO folds)
  - NAS composite (0-8, 9 levels, N=660, 5 LOCO folds)
  - Disease state (healthy/NAFL/NASH, 3 levels, N=1444)

Key outputs:
  - multiprogram_results.csv: per-fold, per-target, per-modality metrics
  - multiprogram_ablation.csv: modality ablation summary
  - multiprogram_random_baselines.csv: per-target random gene baselines
  - multiprogram_concept_attribution.csv: concept → target attribution matrix
  - multiprogram_predictions.csv: per-sample predicted scores

SLURM:
  sbatch --job-name=190_multiprogram \\
         --partition=cpu --cpus-per-task=16 --mem=120G --time=48:00:00 \\
         --output=logs/190_%j.out --error=logs/190_%j.err \\
         --wrap="bash -c 'eval \"\\$(micromamba shell hook --shell bash)\" && \\
                 micromamba activate spatial && \\
                 cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts && \\
                 python 190_multiprogram_predictor.py'"
"""

import os
import sys
import time
import json
import warnings
import logging
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import h5py
from scipy import stats

from sklearn.linear_model import ElasticNet, LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import cohen_kappa_score, mean_absolute_error, accuracy_score

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader

from coral_pytorch.layers import CoralLayer
from coral_pytorch.losses import coral_loss as coral_loss_fn

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

SEED = 42
np.random.seed(SEED)
torch.manual_seed(SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(SEED)

# ── Paths ─────────────────────────────────────────────────────────────────
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
INT = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR = os.path.join(INT, "results")
STAGING = os.path.join(RDIR, "staging_classifier")
OUTDIR = os.path.join(RDIR, "multiprogram")
os.makedirs(OUTDIR, exist_ok=True)

H5_PATH = os.path.join(STAGING, "prepared_data.h5")
META_PATH = os.path.join(STAGING, "modeling_metadata.csv")
DECONV_PATH = os.path.join(RDIR, "deconvolution/bayesprism/unified_bayesprism_proportions.csv")
SSGSEA_PATH = os.path.join(STAGING, "pathway_scores_ssgsea.csv")
TF_PATH = os.path.join(STAGING, "tf_activity_features.csv")
COLOC_PATH = os.path.join(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
TRANSITION_PATH = os.path.join(OUTDIR, "pseudotime_transition_scores.csv")
# Full-universe z-scored expression (~34k genes), the source for the
# leakage-free random-gene null (excludes the supervised-screened 3k pool).
FULL_EXPR_RDS = os.path.join(STAGING, "zscore_expression_full.rds")
FULL_EXPR_CSV = os.path.join(OUTDIR, "_zscore_expression_full.csv")

RSCRIPT = "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
MISSING = -1

# ── Logging ───────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(os.path.join(OUTDIR, "190_multiprogram.log")),
    ],
)
log = logging.getLogger(__name__)

print("=" * 70)
print("190: Multi-Program Decomposition Predictor")
print("=" * 70)
print(f"Device : {device}")
print(f"Output : {OUTDIR}")
print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print()

# ── Task configuration ────────────────────────────────────────────────────
# target_name -> (n_classes, is_ordinal, label_col, fold_col, datasets)
NAS_DATASETS = ["GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066"]
FIB_DATASETS = NAS_DATASETS + ["GSE240729"]

TASK_CFG = {
    "fibrosis": (5, True, "fib_stage", "loco_fold_fibrosis", FIB_DATASETS),
    "nas_composite": (9, True, "nas_score", "loco_fold_nas", NAS_DATASETS),
    "severity": (4, True, "severity4", "loco_fold_fibrosis", FIB_DATASETS),
}
# severity4: 0=Control, 1=NAFL, 2=Borderline, 3=NASH (952 valid, 492 missing)

# ── 1. Load data ──────────────────────────────────────────────────────────
log.info("Loading data...")

# Metadata
meta = pd.read_csv(META_PATH)
sample_ids = meta["sample_id"].values
datasets = meta["dataset"].values
n_samples = len(sample_ids)
log.info(f"  Metadata: {n_samples} samples")

# Expression from H5
with h5py.File(H5_PATH, "r") as h5:
    expr_samples = [s.decode() if isinstance(s, bytes) else s for s in h5["sample_ids"][:]]
    gene_names = [g.decode() if isinstance(g, bytes) else g for g in h5["gene_names"][:]]
    expr_mat = h5["zscore_expression"][:].T  # (genes x samples) -> (samples x genes)
log.info(f"  Expression: {expr_mat.shape[0]} samples x {expr_mat.shape[1]} genes")

# Align expression to metadata sample order
expr_order = {s: i for i, s in enumerate(expr_samples)}
expr_idx = [expr_order[s] for s in sample_ids if s in expr_order]
expr_mat = expr_mat[expr_idx]
log.info(f"  Aligned expression: {expr_mat.shape}")

# Deconvolution proportions — drop non-numeric columns (sample_id, dataset)
deconv_df = pd.read_csv(DECONV_PATH)
deconv_df = deconv_df.set_index("sample_id")
deconv_numeric = deconv_df.select_dtypes(include=[np.number])
deconv = deconv_numeric.reindex(sample_ids).fillna(0).values
deconv_cols = deconv_numeric.columns.tolist()
log.info(f"  Deconvolution: {deconv.shape[1]} cell types")

# TF activity
if os.path.exists(TF_PATH):
    tf_df = pd.read_csv(TF_PATH, index_col=0)
    tf_mat = tf_df.reindex(sample_ids).fillna(0).values
    tf_cols = tf_df.columns.tolist()
    log.info(f"  TF activity: {tf_mat.shape[1]} TFs")
else:
    # Convert from RDS
    log.info("  Converting TF activity from RDS...")
    import subprocess
    rds_path = TF_PATH.replace(".csv", ".rds")
    subprocess.run([RSCRIPT, "-e", f"""
        library(data.table)
        x <- readRDS("{rds_path}")
        fwrite(as.data.table(x, keep.rownames="sample_id"), "{TF_PATH}")
    """], check=True)
    tf_df = pd.read_csv(TF_PATH, index_col=0)
    tf_mat = tf_df.reindex(sample_ids).fillna(0).values
    tf_cols = tf_df.columns.tolist()
    log.info(f"  TF activity: {tf_mat.shape[1]} TFs")

# ssGSEA pathway scores (transposed CSV: samples x pathways, sample_id as last column)
SSGSEA_TRANS = SSGSEA_PATH.replace(".csv", "_transposed.csv")
if os.path.exists(SSGSEA_TRANS):
    ssgsea_df = pd.read_csv(SSGSEA_TRANS)
    # Drop sample_id column if present, keep only numeric pathway columns
    ssgsea_cols = [c for c in ssgsea_df.columns if c != "sample_id"]
    ssgsea = ssgsea_df[ssgsea_cols].values.astype(np.float32)
else:
    # Fall back to original (pathways x samples) and transpose
    ssgsea_raw = pd.read_csv(SSGSEA_PATH, index_col=None)
    ssgsea = ssgsea_raw.values.T.astype(np.float32)  # samples x pathways
    ssgsea_cols = [f"pathway_{i}" for i in range(ssgsea.shape[1])]
log.info(f"  ssGSEA pathways: {ssgsea.shape[1]} pathways")

# COLOC genes — match via Ensembl base ID (strip version suffix)
coloc_df = pd.read_csv(COLOC_PATH)
coloc_ensembl = coloc_df[coloc_df["coloc_best_pp4"] > 0.5]["ensembl"].dropna().unique().tolist()
# Expression matrix gene names have version suffix (ENSG00000121410.14)
# Strip version for matching
gene_base = {g.split(".")[0]: i for i, g in enumerate(gene_names)}
coloc_gene_idx = [gene_base[e] for e in coloc_ensembl if e in gene_base]
coloc_expr = expr_mat[:, coloc_gene_idx]
log.info(f"  COLOC genes: {len(coloc_ensembl)} Ensembl IDs (PP4>0.5), {coloc_expr.shape[1]} in expression matrix")

# ── Full-universe expression for the leakage-free random-gene null ─────────
# NULL FIX (mega-review A7): `expr_mat` (from prepared_data.h5) is NOT the full
# transcriptome — it is the 3,000 genes that Script 61 pre-selected by max |t|
# across the disease/NAS/fibrosis DE contrasts (a SUPERVISED, label-aware
# screen; Script 61 itself flags it "leaks test-fold info"). Drawing a "random"
# baseline from that pool samples genes already enriched for disease signal, so
# the null overstates chance performance and understates how much the
# supervised model beats random. The corrected null draws from the FULL
# ~34k-gene universe EXCLUDING the supervised-screened 3k pool.
full_expr_mat = None
random_pool_idx = None      # indices into full_expr_mat columns, screened pool removed
full_gene_names = None
if os.path.exists(FULL_EXPR_RDS):
    if not os.path.exists(FULL_EXPR_CSV):
        log.info("  Converting full-universe expression RDS -> CSV (one-time)...")
        import subprocess
        subprocess.run([RSCRIPT, "-e", f"""
            suppressMessages(library(data.table))
            m <- readRDS("{FULL_EXPR_RDS}")            # genes x samples
            dt <- as.data.table(m, keep.rownames="gene")
            fwrite(dt, "{FULL_EXPR_CSV}")
        """], check=True)
    full_df = pd.read_csv(FULL_EXPR_CSV).set_index("gene")  # genes x samples
    # Align columns (samples) to metadata sample order; fill missing with 0
    full_df = full_df.reindex(columns=sample_ids).fillna(0.0)
    full_gene_names = full_df.index.tolist()
    full_expr_mat = full_df.values.T.astype(np.float32)      # samples x genes
    # Supervised-screened pool = the 3k genes used by the real model (`gene_names`),
    # matched on base Ensembl ID (strip version suffix on both sides).
    screened_base = {g.split(".")[0] for g in gene_names}
    full_base = np.array([g.split(".")[0] for g in full_gene_names])
    random_pool_idx = np.where(~np.isin(full_base, list(screened_base)))[0]
    log.info(
        f"  Full universe: {full_expr_mat.shape[1]} genes; "
        f"screened pool excluded = {full_expr_mat.shape[1] - len(random_pool_idx)}; "
        f"random-null pool = {len(random_pool_idx)} genes"
    )
else:
    log.warning(
        f"  Full-universe expression not found at {FULL_EXPR_RDS}; "
        f"random-gene null will FALL BACK to the screened 3k pool (NOT leakage-free)."
    )

# Pseudotime transition scores (from Script 191)
transition_scores = None
if os.path.exists(TRANSITION_PATH):
    trans_df = pd.read_csv(TRANSITION_PATH, index_col=0)
    transition_scores = trans_df.reindex(sample_ids).fillna(0).values
    transition_cols = trans_df.columns.tolist()
    log.info(f"  Transition signatures: {transition_scores.shape[1]} scores")
else:
    log.warning("  Transition signature scores not found — skipping pseudotime features")

# Clinical features
clinical = meta[["age", "sex"]].copy()
clinical["sex"] = clinical["sex"].map({"M": 0, "F": 1}).fillna(0.5)
clinical["age"] = clinical["age"].fillna(clinical["age"].median())
clinical_mat = clinical.values.astype(np.float32)
log.info(f"  Clinical: {clinical_mat.shape[1]} features")

# Labels
labels = {}
for task, (n_cls, is_ord, col, fold_col, ds_list) in TASK_CFG.items():
    y = meta[col].fillna(MISSING).astype(int).values
    labels[task] = y
    n_valid = (y != MISSING).sum()
    log.info(f"  Target '{task}': {n_valid} valid samples, {n_cls} classes")


# ── 2. Feature group definitions ──────────────────────────────────────────
def build_feature_groups():
    """Return dict of modality_name -> (feature_matrix, feature_names)."""
    groups = {
        "expression": (expr_mat, [f"gene_{g}" for g in gene_names]),
        "celltype": (deconv, [f"ct_{c}" for c in deconv_cols]),
        "tf_activity": (tf_mat, [f"tf_{c}" for c in tf_cols]),
        "coloc_genetics": (coloc_expr, [f"coloc_{gene_names[i]}" for i in coloc_gene_idx]),
        "ssgsea": (ssgsea, [f"ssgsea_{c}" for c in ssgsea_cols]),
        "clinical": (clinical_mat, ["age", "sex"]),
    }
    if transition_scores is not None:
        groups["pseudotime"] = (transition_scores, [f"trans_{c}" for c in transition_cols])
    return groups


# ── 3. CORAL ordinal utilities (from Script 106) ─────────────────────────
def labels_to_levels(y, n_classes):
    """Convert integer labels -> binary level indicators for CORAL."""
    batch_size = len(y)
    levels = torch.zeros(batch_size, n_classes - 1, dtype=torch.float32)
    for k in range(n_classes - 1):
        levels[:, k] = (y > k).float()
    return levels


def coral_ordinal_probs(logits, n_classes):
    """Convert CORAL logits -> class probabilities."""
    cumprobs = torch.sigmoid(logits)
    batch = cumprobs.shape[0]
    probs = torch.zeros(batch, n_classes, device=logits.device)
    probs[:, 0] = 1.0 - cumprobs[:, 0]
    for k in range(1, n_classes - 1):
        probs[:, k] = cumprobs[:, k - 1] - cumprobs[:, k]
    probs[:, n_classes - 1] = cumprobs[:, n_classes - 2]
    probs = probs.clamp(min=1e-7)
    probs = probs / probs.sum(dim=1, keepdim=True)
    return probs


# ── 4. LOCO-CV with elastic net (per-modality) ───────────────────────────
def run_elastic_net_loco(X, task_name, modality_name, feature_names=None):
    """Run LOCO-CV elastic net for one target + one modality.

    Returns list of dicts with per-fold metrics.
    """
    n_cls, is_ord, col, fold_col, ds_list = TASK_CFG[task_name]
    y = labels[task_name]
    fold_assignments = meta[fold_col].values

    results = []
    all_preds = []

    for fold_ds in ds_list:
        test_mask = (datasets == fold_ds) & (y != MISSING)
        train_mask = (datasets != fold_ds) & (y != MISSING)
        # Only include samples from relevant datasets
        ds_mask = np.isin(datasets, ds_list)
        test_mask = test_mask & ds_mask
        train_mask = train_mask & ds_mask

        if test_mask.sum() == 0 or train_mask.sum() == 0:
            continue

        X_train, X_test = X[train_mask], X[test_mask]
        y_train, y_test = y[train_mask], y[test_mask]

        # Within-fold feature selection: top 500 by variance
        if X_train.shape[1] > 500:
            var = X_train.var(axis=0)
            top_idx = np.argsort(var)[-500:]
            X_train = X_train[:, top_idx]
            X_test = X_test[:, top_idx]

        # Standardize
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        # Ordinal regression via cumulative logit (sklearn approximation)
        # Use OrdinalEncoder approach: fit K-1 binary classifiers
        if is_ord and n_cls > 2:
            preds = np.zeros(len(y_test))
            for threshold in range(1, n_cls):
                y_bin_train = (y_train >= threshold).astype(int)
                y_bin_test = (y_test >= threshold).astype(int)
                if len(np.unique(y_bin_train)) < 2:
                    continue
                lr = LogisticRegression(
                    penalty="elasticnet", solver="saga", l1_ratio=0.5,
                    max_iter=2000, C=1.0, random_state=SEED
                )
                lr.fit(X_train, y_bin_train)
                preds += lr.predict_proba(X_test)[:, 1]
            preds = np.round(preds).astype(int)
            preds = np.clip(preds, 0, n_cls - 1)
        else:
            lr = LogisticRegression(
                penalty="elasticnet", solver="saga", l1_ratio=0.5,
                max_iter=2000, C=1.0, random_state=SEED, multi_class="multinomial"
            )
            lr.fit(X_train, y_train)
            preds = lr.predict(X_test)

        # Metrics
        qwk = cohen_kappa_score(y_test, preds, weights="quadratic")
        mae = mean_absolute_error(y_test, preds)
        acc = accuracy_score(y_test, preds)

        results.append({
            "task": task_name, "modality": modality_name, "fold": fold_ds,
            "n_train": train_mask.sum(), "n_test": test_mask.sum(),
            "n_features": X_train.shape[1],
            "qwk": qwk, "mae": mae, "accuracy": acc,
        })

        for i, idx in enumerate(np.where(test_mask)[0]):
            all_preds.append({
                "sample_id": sample_ids[idx], "task": task_name,
                "modality": modality_name, "fold": fold_ds,
                "y_true": int(y_test[i]), "y_pred": int(preds[i]),
            })

    return results, all_preds


# ── 5. Random gene baseline ──────────────────────────────────────────────
def run_random_gene_baseline(n_draws=100):
    """Run random 500-gene baseline for each target.

    NULL FIX (mega-review A7): draw the 500 random genes from the FULL
    ~34k-gene universe EXCLUDING the supervised-screened 3k pool, so the
    baseline is a genuine null rather than a draw from the disease-enriched
    pre-selected genes. Falls back to the screened pool only if the full
    universe matrix is unavailable (logged loudly upstream).
    """
    rng = np.random.RandomState(SEED)
    baseline_results = []

    if full_expr_mat is not None and random_pool_idx is not None and len(random_pool_idx) >= 500:
        source_mat = full_expr_mat
        pool = random_pool_idx
        log.info(
            f"Running random gene baselines (100 draws) from leakage-free pool "
            f"({len(pool)} genes, supervised 3k pool excluded)..."
        )
    else:
        source_mat = expr_mat
        pool = np.arange(expr_mat.shape[1])
        log.warning(
            "Running random gene baselines (100 draws) from the SCREENED 3k pool "
            "(full universe unavailable) — this null is NOT leakage-free."
        )

    for draw in range(n_draws):
        gene_idx = rng.choice(pool, size=500, replace=False)
        X_random = source_mat[:, gene_idx]

        for task_name in TASK_CFG:
            results, _ = run_elastic_net_loco(X_random, task_name, f"random_{draw}")
            mean_qwk = np.mean([r["qwk"] for r in results]) if results else np.nan
            baseline_results.append({
                "task": task_name, "draw": draw, "mean_qwk": mean_qwk,
            })

        if (draw + 1) % 10 == 0:
            log.info(f"  Random baseline draw {draw + 1}/{n_draws}")

    return pd.DataFrame(baseline_results)


# ── 6. Concept bottleneck model (multi-task CORAL) ───────────────────────

# Concept group sizes (~50 total concepts)
_CONCEPT_GROUPS = {
    "hallmark_pathway":    14,   # supervised by ssGSEA (top 14 by variance)
    "celltype_proportion":  8,   # supervised by deconvolution proportions
    "tf_activity":         14,   # supervised by TF activity (top 14 by variance)
    "masld_program":        8,   # unsupervised MASLD-specific programs
    "other_biology":        6,   # unsupervised catch-all latent concepts
}
N_CONCEPTS = sum(_CONCEPT_GROUPS.values())  # 50


def build_concept_targets_190(sample_mask, ssgsea_mat, ssgsea_names,
                               deconv_mat, deconv_names, tf_mat_local, tf_names):
    """Build concept supervision targets for a subset of samples.

    Parameters
    ----------
    sample_mask : np.ndarray of bool, shape (n_samples,)
        Mask selecting samples to include (train or all).
    ssgsea_mat, deconv_mat, tf_mat_local : np.ndarray (n_samples, n_features)
        Feature matrices aligned to full sample list.
    *_names : list of str

    Returns
    -------
    concept_targets : np.ndarray shape (mask.sum(), N_CONCEPTS)
        Values in [0, 1] for supervised concepts; -1.0 for unsupervised.
    concept_col_names : list of str
    """
    n = int(sample_mask.sum())
    targets = []
    col_names = []

    # Hallmark pathways (top 14 by variance across masked samples only)
    n_hall = _CONCEPT_GROUPS["hallmark_pathway"]
    if ssgsea_mat.shape[1] > 0:
        var_order = np.argsort(ssgsea_mat[sample_mask].var(axis=0))[::-1][:n_hall]
        for rank, idx in enumerate(var_order):
            col = ssgsea_names[idx] if idx < len(ssgsea_names) else f"ssgsea_{idx}"
            vals = ssgsea_mat[sample_mask, idx].astype(np.float32)
            vmin, vmax = vals.min(), vals.max()
            vals = (vals - vmin) / (vmax - vmin) if vmax > vmin else np.full(n, 0.5, np.float32)
            targets.append(vals)
            col_names.append(f"hallmark_{rank}_{col[:25]}")
    else:
        for i in range(n_hall):
            targets.append(np.full(n, -1.0, np.float32))
            col_names.append(f"hallmark_{i}_unknown")

    # Cell-type proportions (up to 8)
    n_ct = _CONCEPT_GROUPS["celltype_proportion"]
    for i in range(n_ct):
        if i < deconv_mat.shape[1]:
            vals = np.clip(deconv_mat[sample_mask, i], 0.0, 1.0).astype(np.float32)
            name = deconv_names[i] if i < len(deconv_names) else f"ct_{i}"
            targets.append(vals)
            col_names.append(f"celltype_{i}_{name}")
        else:
            targets.append(np.full(n, -1.0, np.float32))
            col_names.append(f"celltype_{i}_pad")

    # TF activity (top 14 by variance across masked samples only)
    n_tf = _CONCEPT_GROUPS["tf_activity"]
    if tf_mat_local.shape[1] > 0:
        var_order = np.argsort(tf_mat_local[sample_mask].var(axis=0))[::-1][:n_tf]
        for rank, idx in enumerate(var_order):
            col = tf_names[idx] if idx < len(tf_names) else f"tf_{idx}"
            vals = tf_mat_local[sample_mask, idx].astype(np.float32)
            vmin, vmax = vals.min(), vals.max()
            vals = (vals - vmin) / (vmax - vmin) if vmax > vmin else np.full(n, 0.5, np.float32)
            targets.append(vals)
            col_names.append(f"tf_{rank}_{col[:25]}")
    else:
        for i in range(n_tf):
            targets.append(np.full(n, -1.0, np.float32))
            col_names.append(f"tf_{i}_unknown")

    # Unsupervised MASLD programs + other (filled with -1 = no supervision)
    for i in range(_CONCEPT_GROUPS["masld_program"]):
        targets.append(np.full(n, -1.0, np.float32))
        col_names.append(f"masld_prog_{i}")
    for i in range(_CONCEPT_GROUPS["other_biology"]):
        targets.append(np.full(n, -1.0, np.float32))
        col_names.append(f"other_{i}")

    return np.column_stack(targets).astype(np.float32), col_names


class _ConceptEncoder190(nn.Module):
    """Input → 256 → 128 → n_concepts with sigmoid."""

    def __init__(self, n_input, n_concepts, hidden_dim=256):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_input, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.BatchNorm1d(hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(hidden_dim // 2, n_concepts),
        )

    def forward(self, x):
        return torch.sigmoid(self.net(x))


class _PredictionLayer190(nn.Module):
    """Multi-task CORAL heads for Script 190's 3 tasks.

    Tasks: fibrosis (5-class), nas_composite (9-class), severity (4-class).
    """

    def __init__(self, n_concepts):
        super().__init__()
        hidden = max(64, n_concepts)
        self.shared = nn.Sequential(
            nn.Linear(n_concepts, hidden),
            nn.BatchNorm1d(hidden),
            nn.ReLU(),
            nn.Dropout(0.2),
        )
        self.head_fibrosis      = CoralLayer(hidden, 5)
        self.head_nas_composite = CoralLayer(hidden, 9)
        self.head_severity      = CoralLayer(hidden, 4)

    def forward(self, concepts):
        z = self.shared(concepts)
        return {
            "fibrosis":      self.head_fibrosis(z),
            "nas_composite": self.head_nas_composite(z),
            "severity":      self.head_severity(z),
        }


class ConceptBottleneckModel190(nn.Module):
    """Concept bottleneck for Script 190: encoder → concepts → task heads."""

    def __init__(self, n_input, n_concepts=N_CONCEPTS, hidden_dim=256):
        super().__init__()
        self.concept_encoder = _ConceptEncoder190(n_input, n_concepts, hidden_dim)
        self.predictor = _PredictionLayer190(n_concepts)

    def forward(self, x):
        concepts = self.concept_encoder(x)
        return concepts, self.predictor(concepts)


class _KendallWeights190(nn.Module):
    """Learnable per-task log-variance uncertainty weights (Kendall 2018)."""

    def __init__(self, n_tasks):
        super().__init__()
        self.log_vars = nn.Parameter(torch.zeros(n_tasks))

    def forward(self, losses):
        total = torch.tensor(0.0, device=self.log_vars.device)
        for i, loss_i in enumerate(losses):
            precision = torch.exp(-self.log_vars[i])
            total = total + precision * loss_i + self.log_vars[i]
        return total


class _ConceptDataset190(Dataset):
    """Samples with expression features, task labels, and concept targets."""

    def __init__(self, X, label_dict, concept_targets):
        self.X = torch.tensor(X, dtype=torch.float32)
        self.labels = {k: torch.tensor(v, dtype=torch.long) for k, v in label_dict.items()}
        self.concept_targets = torch.tensor(concept_targets, dtype=torch.float32)
        self.concept_mask = (self.concept_targets >= 0).float()

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        return {
            "X": self.X[idx],
            "labels": {k: v[idx] for k, v in self.labels.items()},
            "concept_targets": self.concept_targets[idx],
            "concept_mask": self.concept_mask[idx],
        }


def _masked_coral_loss_190(logits, targets, n_classes):
    """CORAL loss with -1 label masking."""
    mask = targets >= 0
    if mask.sum() == 0:
        return torch.tensor(0.0, device=logits.device, requires_grad=True)
    logits_v = logits[mask]
    targets_v = targets[mask]
    levels = torch.zeros(len(targets_v), n_classes - 1, dtype=torch.float32, device=logits.device)
    for k in range(n_classes - 1):
        levels[:, k] = (targets_v > k).float()
    return coral_loss_fn(logits_v, levels)


def _train_concept_fold(X_train, lab_train, ct_train,
                        X_val, lab_val, ct_val,
                        n_features, fold_name,
                        phase1_epochs=100, phase2_epochs=200, finetune_epochs=100,
                        lr=5e-4, batch_size=128, patience=30):
    """Three-phase concept bottleneck training for one LOCO fold.

    Phase 1: Concept encoder supervised by biological concept targets (MSE).
    Phase 2: Prediction heads trained with frozen concept encoder (CORAL).
    Phase 3: End-to-end fine-tune at 0.1x LR.

    Returns
    -------
    predictions : dict task -> np.ndarray (n_val, n_classes) probabilities
    val_concepts : np.ndarray (n_val, N_CONCEPTS) concept activations
    model : trained ConceptBottleneckModel190
    """
    train_ds = _ConceptDataset190(X_train, lab_train, ct_train)
    val_ds   = _ConceptDataset190(X_val,   lab_val,   ct_val)
    train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True,
                              drop_last=False, num_workers=0)
    val_loader   = DataLoader(val_ds,   batch_size=batch_size, shuffle=False,
                              num_workers=0)

    model = ConceptBottleneckModel190(
        n_input=n_features,
        n_concepts=N_CONCEPTS,
        hidden_dim=min(512, max(256, n_features // 2)),
    ).to(device)

    # 3 tasks in Script 190
    _TASKS_190 = [("fibrosis", 5), ("nas_composite", 9), ("severity", 4)]
    uncertainty = _KendallWeights190(len(_TASKS_190)).to(device)

    # ---- Phase 1: concept supervision ----
    opt1 = optim.Adam(model.concept_encoder.parameters(), lr=lr, weight_decay=1e-4)
    sch1 = optim.lr_scheduler.CosineAnnealingLR(opt1, T_max=phase1_epochs, eta_min=1e-6)
    best_c_loss, best_c_state, p1_pat = float("inf"), None, 0

    for epoch in range(phase1_epochs):
        model.concept_encoder.train()
        for batch in train_loader:
            X_b = batch["X"].to(device)
            ct  = batch["concept_targets"].to(device)
            cm  = batch["concept_mask"].to(device)
            concepts = model.concept_encoder(X_b)
            diff = ((concepts - ct) ** 2) * cm
            loss = diff.sum() / cm.sum().clamp(min=1)
            opt1.zero_grad(); loss.backward()
            nn.utils.clip_grad_norm_(model.concept_encoder.parameters(), 5.0)
            opt1.step()
        sch1.step()

        model.concept_encoder.eval()
        val_cl = 0.0
        with torch.no_grad():
            for batch in val_loader:
                X_b = batch["X"].to(device)
                ct  = batch["concept_targets"].to(device)
                cm  = batch["concept_mask"].to(device)
                concepts = model.concept_encoder(X_b)
                diff = ((concepts - ct) ** 2) * cm
                val_cl += (diff.sum() / cm.sum().clamp(min=1)).item()
        val_cl /= max(len(val_loader), 1)

        if val_cl < best_c_loss:
            best_c_loss  = val_cl
            best_c_state = {k: v.cpu().clone() for k, v in model.concept_encoder.state_dict().items()}
            p1_pat = 0
        else:
            p1_pat += 1
        if p1_pat >= patience:
            break
        if (epoch + 1) % 25 == 0 or epoch == 0:
            log.info(f"    [{fold_name}] P1 e={epoch+1}: val_concept={val_cl:.4f}")

    if best_c_state:
        model.concept_encoder.load_state_dict(best_c_state)
    model.concept_encoder.to(device)

    # ---- Phase 2: prediction heads (frozen encoder) ----
    for p in model.concept_encoder.parameters():
        p.requires_grad = False
    opt2_params = list(model.predictor.parameters()) + list(uncertainty.parameters())
    opt2 = optim.Adam(opt2_params, lr=lr, weight_decay=1e-4)
    sch2 = optim.lr_scheduler.CosineAnnealingLR(opt2, T_max=phase2_epochs, eta_min=1e-6)
    best_v2, best_p2_state, best_unc_state, p2_pat = float("inf"), None, None, 0

    for epoch in range(phase2_epochs):
        model.predictor.train()
        model.concept_encoder.eval()
        for batch in train_loader:
            X_b = batch["X"].to(device)
            y_b = {k: v.to(device) for k, v in batch["labels"].items()}
            with torch.no_grad():
                concepts = model.concept_encoder(X_b)
            task_logits = model.predictor(concepts)
            task_losses = [
                _masked_coral_loss_190(task_logits[t], y_b[t], nc)
                for t, nc in _TASKS_190
            ]
            loss = uncertainty(task_losses)
            opt2.zero_grad(); loss.backward()
            nn.utils.clip_grad_norm_(opt2_params, 5.0)
            opt2.step()
        sch2.step()

        model.eval()
        val_total = 0.0
        with torch.no_grad():
            for batch in val_loader:
                X_b = batch["X"].to(device)
                y_b = {k: v.to(device) for k, v in batch["labels"].items()}
                concepts = model.concept_encoder(X_b)
                task_logits = model.predictor(concepts)
                for t, nc in _TASKS_190:
                    val_total += _masked_coral_loss_190(task_logits[t], y_b[t], nc).item()
        val_total /= max(len(val_loader), 1)

        if val_total < best_v2:
            best_v2       = val_total
            best_p2_state = {k: v.cpu().clone() for k, v in model.predictor.state_dict().items()}
            best_unc_state = {k: v.cpu().clone() for k, v in uncertainty.state_dict().items()}
            p2_pat = 0
        else:
            p2_pat += 1
        if p2_pat >= patience:
            break
        if (epoch + 1) % 50 == 0 or epoch == 0:
            log.info(f"    [{fold_name}] P2 e={epoch+1}: val={val_total:.4f}")

    if best_p2_state:
        model.predictor.load_state_dict(best_p2_state)
    if best_unc_state:
        uncertainty.load_state_dict(best_unc_state)
    model.predictor.to(device)
    uncertainty.to(device)

    # ---- Phase 3: end-to-end fine-tune ----
    for p in model.concept_encoder.parameters():
        p.requires_grad = True
    all_params = list(model.parameters()) + list(uncertainty.parameters())
    opt3 = optim.Adam(all_params, lr=lr * 0.1, weight_decay=1e-4)
    sch3 = optim.lr_scheduler.CosineAnnealingLR(opt3, T_max=finetune_epochs, eta_min=1e-7)
    best_v3, best_ft_state, p3_pat = float("inf"), None, 0

    for epoch in range(finetune_epochs):
        model.train()
        for batch in train_loader:
            X_b = batch["X"].to(device)
            y_b = {k: v.to(device) for k, v in batch["labels"].items()}
            ct  = batch["concept_targets"].to(device)
            cm  = batch["concept_mask"].to(device)
            concepts, task_logits = model(X_b)
            task_losses = [
                _masked_coral_loss_190(task_logits[t], y_b[t], nc)
                for t, nc in _TASKS_190
            ]
            loss = uncertainty(task_losses)
            # Concept alignment regularisation (weight=0.1)
            diff = ((concepts - ct) ** 2) * cm
            if cm.sum() > 0:
                loss = loss + 0.1 * diff.sum() / cm.sum()
            opt3.zero_grad(); loss.backward()
            nn.utils.clip_grad_norm_(all_params, 5.0)
            opt3.step()
        sch3.step()

        model.eval()
        val_total = 0.0
        with torch.no_grad():
            for batch in val_loader:
                X_b = batch["X"].to(device)
                y_b = {k: v.to(device) for k, v in batch["labels"].items()}
                concepts, task_logits = model(X_b)
                for t, nc in _TASKS_190:
                    val_total += _masked_coral_loss_190(task_logits[t], y_b[t], nc).item()
        val_total /= max(len(val_loader), 1)

        if val_total < best_v3:
            best_v3      = val_total
            best_ft_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            p3_pat = 0
        else:
            p3_pat += 1
        if p3_pat >= patience:
            break
        if (epoch + 1) % 25 == 0 or epoch == 0:
            log.info(f"    [{fold_name}] P3 e={epoch+1}: val={val_total:.4f}")

    if best_ft_state:
        model.load_state_dict(best_ft_state)
    model.to(device)

    # ---- Predict on val set ----
    model.eval()
    all_logits = {"fibrosis": [], "nas_composite": [], "severity": []}
    all_concept_acts = []

    with torch.no_grad():
        for batch in val_loader:
            X_b = batch["X"].to(device)
            concepts, task_logits = model(X_b)
            all_concept_acts.append(concepts.cpu().numpy())
            for t in all_logits:
                all_logits[t].append(task_logits[t].cpu())

    predictions = {}
    task_n_classes = {"fibrosis": 5, "nas_composite": 9, "severity": 4}
    for t, nc in task_n_classes.items():
        cat = torch.cat(all_logits[t], dim=0)
        cumprobs = torch.sigmoid(cat)
        batch_sz  = cumprobs.shape[0]
        probs = torch.zeros(batch_sz, nc)
        probs[:, 0] = 1.0 - cumprobs[:, 0]
        for k in range(1, nc - 1):
            probs[:, k] = cumprobs[:, k - 1] - cumprobs[:, k]
        probs[:, nc - 1] = cumprobs[:, nc - 2]
        probs = probs.clamp(min=1e-7)
        probs = probs / probs.sum(dim=1, keepdim=True)
        predictions[t] = probs.numpy()

    val_concepts = np.vstack(all_concept_acts)
    return predictions, val_concepts, model


def _compute_concept_attribution(model, X_train, concept_col_names):
    """Compute concept→target attribution via gradient of CORAL logit mean.

    For each (concept_i, target_j) pair: mean |d(sum logit_j) / d(concept_i)|
    averaged over a random subsample of training samples.

    Returns
    -------
    attr : np.ndarray (N_CONCEPTS, n_tasks)
    """
    _TASKS_ATT = [("fibrosis", 5), ("nas_composite", 9), ("severity", 4)]
    model.eval()
    n_sub = min(256, len(X_train))
    rng = np.random.RandomState(SEED)
    idx = rng.choice(len(X_train), n_sub, replace=False)
    X_sub = torch.tensor(X_train[idx], dtype=torch.float32).to(device)

    attr_accum = np.zeros((N_CONCEPTS, len(_TASKS_ATT)), dtype=np.float32)

    for b_start in range(0, n_sub, 64):
        X_b = X_sub[b_start: b_start + 64]
        concepts = model.concept_encoder(X_b)  # (batch, N_CONCEPTS)
        concepts = concepts.detach().requires_grad_(True)

        task_logits = model.predictor(concepts)

        for j, (t, _nc) in enumerate(_TASKS_ATT):
            # Use sum of all CORAL logits as scalar signal
            scalar = task_logits[t].sum()
            grad = torch.autograd.grad(
                scalar, concepts, retain_graph=(j < len(_TASKS_ATT) - 1)
            )[0]  # (batch, N_CONCEPTS)
            attr_accum[:, j] += grad.abs().mean(dim=0).detach().cpu().numpy()

    attr_accum /= max(1, (n_sub + 63) // 64)
    return attr_accum


def run_concept_bottleneck_loco():
    """Run concept bottleneck LOCO-CV for all 3 targets jointly.

    Uses the same feature matrix as elastic net (combined all modalities),
    plus biological concept supervision from ssGSEA, deconv, and TF activity.

    Returns
    -------
    results_df   : pd.DataFrame  per-fold × per-target QWK
    attribution_df : pd.DataFrame  concept × target attribution (averaged across folds)
    """
    log.info("Building combined feature matrix for concept bottleneck...")
    feature_groups = build_feature_groups()
    X_combined = np.hstack([v[0] for v in feature_groups.values()]).astype(np.float32)
    n_features = X_combined.shape[1]
    log.info(f"  Combined features: {n_features}")

    # Concept target components (full-matrix, will be sliced per fold)
    deconv_cols_list = [f"ct_{c}" for c in deconv_cols]

    all_results = []
    all_attr = np.zeros((N_CONCEPTS, 3), dtype=np.float32)
    n_folds_done = 0
    concept_col_names = None  # captured on first fold

    # Iterate over all LOCO folds across all tasks
    # We use the union of fold datasets (FIB_DATASETS covers NAS_DATASETS too)
    # and skip folds where a specific task has no valid labels in the test set.
    for fold_ds in FIB_DATASETS:
        ds_mask_all = np.isin(datasets, FIB_DATASETS)
        test_mask_base = (datasets == fold_ds) & ds_mask_all
        train_mask_base = (datasets != fold_ds) & ds_mask_all

        if test_mask_base.sum() == 0 or train_mask_base.sum() == 0:
            continue

        fold_name = f"cbm_{fold_ds}"
        log.info(f"  Fold: {fold_name} (train={train_mask_base.sum()}, test={test_mask_base.sum()})")

        X_train = X_combined[train_mask_base]
        X_test  = X_combined[test_mask_base]

        # Standardize (fit on train)
        scaler = StandardScaler()
        X_train_sc = scaler.fit_transform(X_train)
        X_test_sc  = scaler.transform(X_test)

        # Build concept targets
        ct_matrix, concept_col_names = build_concept_targets_190(
            sample_mask=train_mask_base,
            ssgsea_mat=ssgsea, ssgsea_names=ssgsea_cols,
            deconv_mat=deconv, deconv_names=deconv_cols,
            tf_mat_local=tf_mat, tf_names=tf_cols,
        )
        # Validation concept targets (for val_ds creation — values not used in loss)
        ct_val_matrix, _ = build_concept_targets_190(
            sample_mask=test_mask_base,
            ssgsea_mat=ssgsea, ssgsea_names=ssgsea_cols,
            deconv_mat=deconv, deconv_names=deconv_cols,
            tf_mat_local=tf_mat, tf_names=tf_cols,
        )

        # Assemble per-task labels (using full label arrays, masked via sample_mask)
        lab_train = {t: labels[t][train_mask_base] for t in TASK_CFG}
        lab_test  = {t: labels[t][test_mask_base]  for t in TASK_CFG}

        # Train
        predictions, val_concepts, trained_model = _train_concept_fold(
            X_train_sc, lab_train, ct_matrix,
            X_test_sc,  lab_test,  ct_val_matrix,
            n_features=n_features, fold_name=fold_name,
        )

        # Evaluate per task
        task_n_classes = {"fibrosis": 5, "nas_composite": 9, "severity": 4}
        for task_name, n_cls in task_n_classes.items():
            y_true = lab_test[task_name]
            valid  = y_true >= 0
            if valid.sum() < 5:
                continue
            y_v   = y_true[valid].astype(int)
            probs = predictions[task_name][valid]
            y_pred = probs.argmax(axis=1)
            y_pred = np.clip(y_pred, 0, n_cls - 1)
            qwk = cohen_kappa_score(y_v, y_pred, weights="quadratic") if len(np.unique(y_v)) > 1 else np.nan
            mae = float(np.mean(np.abs(y_v - y_pred)))
            acc = float(accuracy_score(y_v, y_pred))
            all_results.append({
                "task": task_name, "fold": fold_ds,
                "n_train": int(train_mask_base.sum()),
                "n_test": int(valid.sum()),
                "qwk": qwk, "mae": mae, "accuracy": acc,
            })
            log.info(f"    {task_name}: QWK={qwk:.3f} MAE={mae:.2f} Acc={acc:.3f}")

        # Attribution
        fold_attr = _compute_concept_attribution(trained_model, X_train_sc, concept_col_names)
        all_attr += fold_attr
        n_folds_done += 1

    results_df = pd.DataFrame(all_results)

    # Average attribution across folds
    if n_folds_done > 0:
        all_attr /= n_folds_done
    task_names_attr = ["fibrosis", "nas_composite", "severity"]
    attribution_df = pd.DataFrame(
        all_attr,
        index=concept_col_names if concept_col_names else [f"concept_{i}" for i in range(N_CONCEPTS)],
        columns=task_names_attr,
    )
    attribution_df.index.name = "concept"

    return results_df, attribution_df


# ── 6. Modality ablation ─────────────────────────────────────────────────
def run_modality_ablation():
    """Run elastic net LOCO-CV for each modality × each target."""
    log.info("Running modality ablation...")
    feature_groups = build_feature_groups()

    all_results = []
    all_preds = []

    for mod_name, (X_mod, feat_names) in feature_groups.items():
        log.info(f"  Modality: {mod_name} ({X_mod.shape[1]} features)")
        for task_name in TASK_CFG:
            results, preds = run_elastic_net_loco(X_mod, task_name, mod_name, feat_names)
            all_results.extend(results)
            all_preds.extend(preds)

    # Combined (all features)
    log.info("  Modality: combined (all features)")
    combined_parts = [v[0] for v in feature_groups.values()]
    X_combined = np.hstack(combined_parts)
    for task_name in TASK_CFG:
        results, preds = run_elastic_net_loco(X_combined, task_name, "combined")
        all_results.extend(results)
        all_preds.extend(preds)

    return pd.DataFrame(all_results), pd.DataFrame(all_preds)


# ── 7. Main execution ────────────────────────────────────────────────────
if __name__ == "__main__":
    # Part A: Modality ablation with elastic net
    log.info("=" * 50)
    log.info("PART A: Modality ablation (elastic net)")
    log.info("=" * 50)
    ablation_results, ablation_preds = run_modality_ablation()
    ablation_results.to_csv(os.path.join(OUTDIR, "multiprogram_ablation_raw.csv"), index=False)
    ablation_preds.to_csv(os.path.join(OUTDIR, "multiprogram_predictions_enet.csv"), index=False)

    # Summarize ablation
    summary = ablation_results.groupby(["task", "modality"]).agg(
        mean_qwk=("qwk", "mean"), std_qwk=("qwk", "std"),
        mean_mae=("mae", "mean"), mean_acc=("accuracy", "mean"),
        n_folds=("qwk", "count"),
    ).reset_index()
    summary.to_csv(os.path.join(OUTDIR, "multiprogram_ablation.csv"), index=False)
    log.info(f"Ablation summary:\n{summary.to_string()}")

    # Part B: Random gene baselines
    log.info("=" * 50)
    log.info("PART B: Random gene baselines")
    log.info("=" * 50)
    random_baselines = run_random_gene_baseline(n_draws=100)
    random_baselines.to_csv(os.path.join(OUTDIR, "multiprogram_random_baselines.csv"), index=False)

    # Compute p-values (expression vs random per target)
    for task_name in TASK_CFG:
        expr_qwk = summary[(summary["task"] == task_name) & (summary["modality"] == "expression")]["mean_qwk"].values
        rand_qwks = random_baselines[random_baselines["task"] == task_name]["mean_qwk"].values
        if len(expr_qwk) > 0 and len(rand_qwks) > 0:
            p_val = (np.sum(rand_qwks >= expr_qwk[0]) + 1) / (len(rand_qwks) + 1)
            log.info(f"  {task_name}: expression QWK={expr_qwk[0]:.3f}, random mean={rand_qwks.mean():.3f}, p={p_val:.3f}")

    # Part C: Concept bottleneck LOCO-CV
    log.info("=" * 50)
    log.info("PART C: Concept bottleneck LOCO-CV (multi-task CORAL)")
    log.info("=" * 50)
    cbm_results, cbm_attribution = run_concept_bottleneck_loco()

    cbm_results.to_csv(os.path.join(OUTDIR, "multiprogram_concept_results.csv"), index=False)
    cbm_attribution.to_csv(os.path.join(OUTDIR, "multiprogram_concept_attribution.csv"))

    if not cbm_results.empty:
        cbm_summary = cbm_results.groupby("task").agg(
            mean_qwk=("qwk", "mean"), std_qwk=("qwk", "std"),
            mean_mae=("mae", "mean"), n_folds=("qwk", "count"),
        ).reset_index()
        log.info(f"Concept bottleneck summary:\n{cbm_summary.to_string()}")

    log.info(f"\nCompleted at {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log.info(f"Results in: {OUTDIR}")
