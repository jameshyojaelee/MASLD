#!/usr/bin/env python3
"""
106_concept_bottleneck.py
COMPASS-style concept bottleneck model for interpretable MASLD staging.

Architecture:
  Input: Multi-scale features (~560 features)
    - ssGSEA pathway scores (~200) from Script 103
    - WGCNA eigengenes (~20-30) from Script 104
    - TF activity scores (~70 significant) from Script 105
    - Deconv proportions (4-8) from BayesPrism
    - VAE embeddings (64) from Script 69
    - Sex (1)
        |
        v
  Concept Layer (Linear -> Sigmoid, ~50 concepts):
    - 14 Hallmark pathway activities (supervised by ssGSEA scores)
    - 8 cell-type proportions (supervised by deconv)
    - 14 key TF activities (supervised by decoupleR)
    - 8 MASLD-specific programs (steatosis/inflammation/fibrosis signatures)
    - 6 other biological features
        |
        v
  Prediction Layer (CORAL ordinal heads, 7 tasks):
    - Disease binary, NAS 3-class, NAS 9-class, Steatosis,
      Inflammation, Ballooning, Fibrosis

Training:
  - Two-phase:
      Phase 1: Train concept layer to reconstruct known biological scores
      Phase 2: Train prediction layer using concepts as input (frozen concept encoder)
    Then fine-tune end-to-end with reduced LR.
  - 5-fold NAS LOCO + 6-fold fibrosis LOCO
  - CORAL ordinal loss for ordinal tasks
  - Masked loss for missing labels (-1)
  - Kendall uncertainty weighting for multi-task balancing

Input (all from results/staging_classifier/):
  - pathway_scores_ssgsea.rds  (Script 103; convert via Rscript)
  - wgcna_eigengenes.rds       (Script 104; convert via Rscript)
  - tf_activity_features.rds   (Script 105; convert via Rscript)
  - modeling_metadata.csv       (labels)
  - embeddings_all_samples.csv  (VAE embeddings from Script 69)
  - unified_bayesprism_proportions.csv (deconv from Script 22)
  - prepared_data.h5            (expression for sample ordering)

Output (all to results/staging_classifier/):
  - concept_bottleneck_results.csv     (per-fold per-task metrics)
  - concept_bottleneck_summary.csv     (aggregated)
  - concept_activations.csv            (per-sample concept scores)
  - concept_bottleneck_model.pt        (best model weights)

SLURM: gpu partition, 1xL40S, 8 CPUs, 64GB RAM, 48h
Env:   micromamba activate rapids_singlecell

Usage:
  sbatch --job-name=stg106_concept \
         --partition=gpu --gres=gpu:1 \
         --cpus-per-task=8 --mem=64G --time=48:00:00 \
         --output=logs/106_concept_%j.out \
         --error=logs/106_concept_%j.err \
         --wrap="bash -c 'eval \"\\$(micromamba shell hook --shell bash)\" && \\
                 micromamba activate rapids_singlecell && \\
                 cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts && \\
                 python 106_concept_bottleneck.py'"
"""

import os
import sys
import time
import subprocess
import warnings
import logging
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import h5py

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader

from coral_pytorch.layers import CoralLayer
from coral_pytorch.losses import coral_loss as coral_loss_fn

from sklearn.metrics import roc_auc_score, accuracy_score

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------
SEED = 42
np.random.seed(SEED)
torch.manual_seed(SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(SEED)

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

# Input files
H5_PATH = os.path.join(OUTDIR, "prepared_data.h5")
META_PATH = os.path.join(OUTDIR, "modeling_metadata.csv")
EMBED_PATH = os.path.join(OUTDIR, "embeddings_all_samples.csv")
DECONV_PATH = os.path.join(
    INT, "results/deconvolution/bayesprism/unified_bayesprism_proportions.csv"
)

# RDS inputs (from upstream scripts 103/104/105)
SSGSEA_RDS = os.path.join(OUTDIR, "pathway_scores_ssgsea.rds")
WGCNA_RDS = os.path.join(OUTDIR, "wgcna_eigengenes.rds")
TF_RDS = os.path.join(OUTDIR, "tf_activity_features.rds")

# Outputs
OUT_RESULTS = os.path.join(OUTDIR, "concept_bottleneck_results.csv")
OUT_SUMMARY = os.path.join(OUTDIR, "concept_bottleneck_summary.csv")
OUT_CONCEPTS = os.path.join(OUTDIR, "concept_activations.csv")
OUT_MODEL = os.path.join(OUTDIR, "concept_bottleneck_model.pt")

# Rscript
RSCRIPT = "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"

MISSING = -1

# Device
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(os.path.join(LOGDIR, "106_concept_bottleneck.log")),
    ],
)
log = logging.getLogger(__name__)

print("=" * 70)
print("106: Concept Bottleneck Model for MASLD Staging")
print("=" * 70)
print(f"Device : {device}")
print(f"Output : {OUTDIR}")
print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print()

# =============================================================================
# Task configuration
# =============================================================================
# Task name -> (n_classes, task_weight, is_ordinal)
TASK_CFG = {
    "disease":       (2, 0.1, False),
    "nas_3class":    (3, 0.3, True),
    "nas_9class":    (9, 0.5, True),
    "steatosis":     (4, 2.0, True),
    "inflammation":  (3, 2.0, True),
    "ballooning":    (3, 2.0, True),
    "fibrosis":      (5, 0.3, True),
}

# LOCO fold definitions
NAS_DATASETS = ["GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066"]
FIB_DATASETS = [
    "GSE130970", "GSE135251", "GSE162694", "GSE174478",
    "GSE193066", "GSE240729",
]

# Concept definitions: name -> (type, n_concepts)
# These define the interpretable bottleneck
CONCEPT_GROUPS = {
    "hallmark_pathway":   14,   # supervised by ssGSEA Hallmark scores
    "celltype_proportion": 8,   # supervised by deconvolution
    "tf_activity":        14,   # supervised by decoupleR TF activity
    "masld_program":       8,   # steatosis, inflammation, fibrosis sigs + more
    "other_biology":       6,   # catch-all latent biological concepts
}
N_CONCEPTS = sum(CONCEPT_GROUPS.values())  # 50

# =============================================================================
# RDS conversion utility
# =============================================================================
def rds_to_csv(rds_path, csv_path, force=False):
    """Convert an RDS file (matrix/data.frame) to CSV via Rscript."""
    if os.path.exists(csv_path) and not force:
        log.info(f"  CSV already exists: {csv_path}")
        return csv_path

    if not os.path.exists(rds_path):
        log.warning(f"  RDS not found: {rds_path}")
        return None

    r_code = f"""
obj <- readRDS("{rds_path}")
if (is.matrix(obj)) {{
    # Check if rows are samples (SRR/GSM/ERR IDs) or features (pathways/genes)
    # If rownames look like pathways/genes (not sample IDs), transpose
    rn <- rownames(obj)
    is_sample_rows <- any(grepl("^SRR|^GSM|^ERR|^PRJNA", rn[1:min(5, length(rn))]))
    if (!is_sample_rows && any(grepl("^SRR|^GSM|^ERR|^PRJNA", colnames(obj)[1:min(5, ncol(obj))]))) {{
        cat("Transposing matrix (features x samples -> samples x features)\\n")
        obj <- t(obj)
    }}
    df <- as.data.frame(obj)
    df$sample_id <- rownames(obj)
}} else if (is.data.frame(obj)) {{
    df <- obj
    if (!"sample_id" %in% colnames(df)) {{
        df$sample_id <- rownames(obj)
    }}
}} else {{
    stop("RDS object is not a matrix or data.frame")
}}
write.csv(df, "{csv_path}", row.names=FALSE)
cat("Wrote", nrow(df), "x", ncol(df), "to {csv_path}\\n")
"""
    r_script = csv_path.replace(".csv", "_convert.R")
    with open(r_script, "w") as f:
        f.write(r_code)

    result = subprocess.run(
        [RSCRIPT, r_script],
        capture_output=True, text=True, timeout=300,
    )
    if result.returncode != 0:
        log.warning(f"  Rscript failed: {result.stderr[:500]}")
        return None

    log.info(f"  {result.stdout.strip()}")
    # Clean up temp R script
    if os.path.exists(r_script):
        os.remove(r_script)
    return csv_path


# =============================================================================
# Load and align all feature matrices
# =============================================================================
def load_features():
    """Load all input feature matrices and align to common sample order."""
    log.info("Loading feature matrices...")

    # 1. Sample IDs and metadata from H5 + CSV
    with h5py.File(H5_PATH, "r") as h5:
        sample_ids = np.array(
            [s.decode() if isinstance(s, bytes) else s for s in h5["sample_ids"][:]]
        )
    n_total = len(sample_ids)
    log.info(f"  Total samples (H5): {n_total}")

    meta_df = pd.read_csv(META_PATH)
    meta_df = meta_df.set_index("sample_id").reindex(sample_ids).reset_index()
    assert len(meta_df) == n_total

    # 2. VAE embeddings (always available)
    log.info("  Loading VAE embeddings...")
    vae_df = pd.read_csv(EMBED_PATH)
    vae_df = vae_df.set_index("sample_id").reindex(sample_ids).fillna(0.0)
    vae_cols = [c for c in vae_df.columns if c.startswith("z")]
    log.info(f"    VAE: {len(vae_cols)} dimensions")

    # 3. Deconvolution proportions
    deconv_df = None
    deconv_cols = []
    if os.path.exists(DECONV_PATH):
        log.info("  Loading deconvolution proportions...")
        raw = pd.read_csv(DECONV_PATH)
        # Use bp_ (BayesPrism) columns
        bp_cols = [c for c in raw.columns if c.startswith("bp_")]
        if bp_cols:
            deconv_df = raw.set_index("sample_id")[bp_cols].reindex(sample_ids).fillna(0.0)
            deconv_cols = bp_cols
            log.info(f"    Deconv: {len(deconv_cols)} cell types ({deconv_cols})")
    if deconv_df is None:
        log.warning("  Deconv proportions not available")
        deconv_df = pd.DataFrame(index=sample_ids)

    # 4. ssGSEA pathway scores (from Script 103)
    ssgsea_df = None
    ssgsea_cols = []
    ssgsea_csv = os.path.join(OUTDIR, "pathway_scores_ssgsea.csv")
    csv_path = rds_to_csv(SSGSEA_RDS, ssgsea_csv)
    if csv_path and os.path.exists(csv_path):
        log.info("  Loading ssGSEA pathway scores...")
        raw = pd.read_csv(csv_path)
        if "sample_id" in raw.columns:
            raw = raw.set_index("sample_id")
        ssgsea_df = raw.reindex(sample_ids).fillna(0.0)
        ssgsea_cols = list(ssgsea_df.columns)
        log.info(f"    ssGSEA: {len(ssgsea_cols)} pathways")
    else:
        log.warning("  ssGSEA scores not available — will use expression fallback")
        ssgsea_df = pd.DataFrame(index=sample_ids)

    # 5. WGCNA eigengenes (from Script 104)
    wgcna_df = None
    wgcna_cols = []
    wgcna_csv = os.path.join(OUTDIR, "wgcna_eigengenes.csv")
    csv_path = rds_to_csv(WGCNA_RDS, wgcna_csv)
    if csv_path and os.path.exists(csv_path):
        log.info("  Loading WGCNA eigengenes...")
        raw = pd.read_csv(csv_path)
        if "sample_id" in raw.columns:
            raw = raw.set_index("sample_id")
        wgcna_df = raw.reindex(sample_ids).fillna(0.0)
        wgcna_cols = list(wgcna_df.columns)
        log.info(f"    WGCNA: {len(wgcna_cols)} eigengenes")
    else:
        log.warning("  WGCNA eigengenes not available — will skip")
        wgcna_df = pd.DataFrame(index=sample_ids)

    # 6. TF activity scores (from Script 105)
    tf_df = None
    tf_cols = []
    tf_csv = os.path.join(OUTDIR, "tf_activity_features.csv")
    csv_path = rds_to_csv(TF_RDS, tf_csv)
    if csv_path and os.path.exists(csv_path):
        log.info("  Loading TF activity scores...")
        raw = pd.read_csv(csv_path)
        if "sample_id" in raw.columns:
            raw = raw.set_index("sample_id")
        tf_df = raw.reindex(sample_ids).fillna(0.0)
        tf_cols = list(tf_df.columns)
        log.info(f"    TF activity: {len(tf_cols)} TFs")
    else:
        log.warning("  TF activity scores not available — will skip")
        tf_df = pd.DataFrame(index=sample_ids)

    # 7. Sex as numeric feature
    sex_arr = np.zeros(n_total, dtype=np.float32)
    if "sex" in meta_df.columns:
        sex_arr = (meta_df["sex"].values == "F").astype(np.float32)
    sex_df = pd.DataFrame({"sex_female": sex_arr}, index=sample_ids)

    # 8. If ssGSEA/WGCNA/TF unavailable, supplement with rank expression
    # from H5 to ensure adequate feature dimensionality
    expr_fallback_df = pd.DataFrame(index=sample_ids)
    if len(ssgsea_cols) == 0 and len(wgcna_cols) == 0 and len(tf_cols) == 0:
        log.info("  Loading rank expression as feature fallback...")
        with h5py.File(H5_PATH, "r") as h5:
            rank_expr = h5["rank_expression"][:].T.astype(np.float32)
            gene_names = [
                g.decode() if isinstance(g, bytes) else g
                for g in h5["gene_names"][:]
            ]
        # Use top ~200 most variable genes as proxy features
        gene_var = np.var(rank_expr, axis=0)
        top_idx = np.argsort(gene_var)[-200:]
        top_genes = [gene_names[i] for i in top_idx]
        expr_fallback_df = pd.DataFrame(
            rank_expr[:, top_idx],
            columns=[f"expr_{g}" for g in top_genes],
            index=sample_ids,
        )
        log.info(f"    Expression fallback: {expr_fallback_df.shape[1]} features")

    # Concatenate all features
    feature_dfs = [
        vae_df[vae_cols] if vae_cols else pd.DataFrame(index=sample_ids),
        deconv_df,
        ssgsea_df,
        wgcna_df,
        tf_df,
        sex_df,
        expr_fallback_df,
    ]
    X_df = pd.concat([df for df in feature_dfs if len(df.columns) > 0], axis=1)
    X = X_df.values.astype(np.float32)

    # Feature group tracking (for concept supervision)
    feature_groups = {
        "vae": vae_cols,
        "deconv": deconv_cols,
        "ssgsea": ssgsea_cols,
        "wgcna": wgcna_cols,
        "tf": tf_cols,
        "sex": ["sex_female"],
        "expr_fallback": list(expr_fallback_df.columns),
    }

    log.info(f"\n  Combined feature matrix: {X.shape[0]} samples x {X.shape[1]} features")
    for group, cols in feature_groups.items():
        if cols:
            log.info(f"    {group}: {len(cols)} features")

    # Build label arrays
    log.info("\nBuilding label arrays...")
    labels = {}

    def safe_int(series, valid_range=None):
        arr = series.fillna(MISSING).values.astype(np.int64)
        if valid_range is not None:
            lo, hi = valid_range
            arr[(arr < lo) | (arr > hi)] = MISSING
        return arr

    labels["disease"] = safe_int(meta_df["is_disease"], (0, 1))
    labels["nas_3class"] = safe_int(meta_df["nas_3class"], (0, 2))
    nas_raw = meta_df["nas_score"].fillna(-1).values.astype(np.int64)
    nas_raw[nas_raw > 8] = 8
    nas_raw[nas_raw < 0] = MISSING
    labels["nas_9class"] = nas_raw
    labels["steatosis"] = safe_int(meta_df["steatosis_grade"], (0, 3))
    labels["inflammation"] = safe_int(meta_df["lobular_inflammation_grade"], (0, 2))
    labels["ballooning"] = safe_int(meta_df["ballooning_grade"], (0, 2))
    labels["fibrosis"] = safe_int(meta_df["fib_stage"], (0, 4))

    for t, arr in labels.items():
        n_valid = int(np.sum(arr >= 0))
        log.info(f"  {t:20s}: {n_valid:5d} labels")

    datasets_arr = meta_df["dataset"].values

    # Build concept supervision targets
    # These are the "ground truth" concept values used in Phase 1 training
    concept_targets = build_concept_targets(
        sample_ids, ssgsea_df, ssgsea_cols,
        deconv_df, deconv_cols,
        tf_df, tf_cols,
    )

    return X, labels, sample_ids, datasets_arr, feature_groups, concept_targets


def build_concept_targets(sample_ids, ssgsea_df, ssgsea_cols,
                          deconv_df, deconv_cols, tf_df, tf_cols):
    """Build RAW (unnormalized) target values for concept supervision.

    Returns dict mapping concept_name -> np.array of shape (n_samples,).
    Hallmark and TF values are raw floats; cell-type proportions are clipped
    to [0, 1]; unsupervised concepts are -1.0 sentinels.

    Normalization to [0, 1] is deferred to normalize_concept_targets() which
    must be called per LOCO fold using train-fold statistics.
    """
    targets = {}
    n = len(sample_ids)

    # Hallmark pathway concepts (top 14 by variance)
    if ssgsea_cols:
        var_order = np.argsort(
            ssgsea_df[ssgsea_cols].var().values
        )[::-1][:CONCEPT_GROUPS["hallmark_pathway"]]
        for i, idx in enumerate(var_order):
            col = ssgsea_cols[idx]
            vals = ssgsea_df[col].values.astype(np.float32)
            targets[f"hallmark_{i}_{col[:30]}"] = vals
    else:
        for i in range(CONCEPT_GROUPS["hallmark_pathway"]):
            targets[f"hallmark_{i}_unknown"] = np.full(n, -1.0, dtype=np.float32)

    # Cell-type proportion concepts (already in [0, 1])
    if deconv_cols:
        for i, col in enumerate(deconv_cols[:CONCEPT_GROUPS["celltype_proportion"]]):
            vals = deconv_df[col].values.astype(np.float32)
            vals = np.clip(vals, 0, 1)
            targets[f"celltype_{i}_{col}"] = vals
        # Pad if fewer than expected
        for i in range(len(deconv_cols), CONCEPT_GROUPS["celltype_proportion"]):
            targets[f"celltype_{i}_pad"] = np.full(n, -1.0, dtype=np.float32)
    else:
        for i in range(CONCEPT_GROUPS["celltype_proportion"]):
            targets[f"celltype_{i}_unknown"] = np.full(n, -1.0, dtype=np.float32)

    # TF activity concepts (top 14 by variance)
    if tf_cols:
        var_order = np.argsort(
            tf_df[tf_cols].var().values
        )[::-1][:CONCEPT_GROUPS["tf_activity"]]
        for i, idx in enumerate(var_order):
            col = tf_cols[idx]
            vals = tf_df[col].values.astype(np.float32)
            targets[f"tf_{i}_{col[:30]}"] = vals
    else:
        for i in range(CONCEPT_GROUPS["tf_activity"]):
            targets[f"tf_{i}_unknown"] = np.full(n, -1.0, dtype=np.float32)

    # MASLD program concepts + other (unsupervised — no targets)
    for i in range(CONCEPT_GROUPS["masld_program"]):
        targets[f"masld_prog_{i}"] = np.full(n, -1.0, dtype=np.float32)
    for i in range(CONCEPT_GROUPS["other_biology"]):
        targets[f"other_{i}"] = np.full(n, -1.0, dtype=np.float32)

    concept_names = list(targets.keys())
    log.info(f"  Concept targets: {len(concept_names)} concepts")
    n_supervised = sum(1 for v in targets.values() if v[0] >= 0)
    log.info(f"  Supervised concepts: {n_supervised}")

    return targets


def normalize_concept_targets(raw_targets, train_mask):
    """Min-max normalize concept targets using train-fold statistics only.

    Parameters
    ----------
    raw_targets : dict[str, np.ndarray]
        Raw concept target arrays from build_concept_targets(), shape (n_samples,).
    train_mask : np.ndarray of bool, shape (n_samples,)
        Boolean mask identifying training samples for this LOCO fold.

    Returns
    -------
    dict[str, np.ndarray]
        Normalized targets (same shape). Hallmark and TF concepts are
        min-max scaled to [0, 1] using train-fold min/max. Cell-type
        proportions and unsupervised sentinels (-1.0) are passed through.
    """
    normed = {}
    n = len(next(iter(raw_targets.values())))
    for key, vals in raw_targets.items():
        # Skip unsupervised concepts (sentinel -1.0) and cell-type proportions
        # (already in [0, 1]). Normalize hallmark_* and tf_* concepts.
        if vals[0] < 0 or key.startswith("celltype_"):
            normed[key] = vals.copy()
            continue
        train_vals = vals[train_mask]
        vmin, vmax = train_vals.min(), train_vals.max()
        if vmax > vmin:
            normed[key] = ((vals - vmin) / (vmax - vmin)).astype(np.float32)
        else:
            normed[key] = np.full(n, 0.5, dtype=np.float32)
    return normed


# =============================================================================
# CORAL ordinal utilities
# =============================================================================
def labels_to_levels(y, n_classes):
    """Convert integer labels -> binary level indicators for CORAL.

    For label y and K classes: levels[k] = 1 if y > k, else 0 (K-1 levels).
    Missing labels (y == -1) -> all zeros (masked out later).
    """
    batch_size = len(y)
    levels = torch.zeros(batch_size, n_classes - 1, dtype=torch.float32)
    for k in range(n_classes - 1):
        levels[:, k] = (y > k).float()
    return levels


def coral_ordinal_probs(logits, n_classes):
    """Convert CORAL logits -> class probabilities.

    P(Y > k) = sigmoid(logit_k)
    P(Y = k) = P(Y > k-1) - P(Y > k)
    """
    cumprobs = torch.sigmoid(logits)  # (batch, K-1)
    batch = cumprobs.shape[0]
    probs = torch.zeros(batch, n_classes, device=logits.device)

    # P(Y = 0) = 1 - P(Y > 0)
    probs[:, 0] = 1.0 - cumprobs[:, 0]

    # P(Y = k) = P(Y > k-1) - P(Y > k) for 1 <= k < K-1
    for k in range(1, n_classes - 1):
        probs[:, k] = cumprobs[:, k - 1] - cumprobs[:, k]

    # P(Y = K-1) = P(Y > K-2)
    probs[:, n_classes - 1] = cumprobs[:, n_classes - 2]

    # Clamp and renormalize for numerical safety
    probs = probs.clamp(min=1e-7)
    probs = probs / probs.sum(dim=1, keepdim=True)
    return probs


def masked_coral_loss(logits, targets, n_classes):
    """CORAL loss with masking for missing labels (-1)."""
    mask = targets >= 0
    if mask.sum() == 0:
        return torch.tensor(0.0, device=logits.device, requires_grad=True)

    logits_v = logits[mask]
    targets_v = targets[mask]
    levels = labels_to_levels(targets_v, n_classes).to(logits.device)
    return coral_loss_fn(logits_v, levels)


def masked_cross_entropy(logits, targets, weight=None):
    """Standard CE with masking for missing labels (-1)."""
    mask = targets >= 0
    if mask.sum() == 0:
        return torch.tensor(0.0, device=logits.device, requires_grad=True)
    return F.cross_entropy(logits[mask], targets[mask], weight=weight)


# =============================================================================
# Model architecture
# =============================================================================
class ConceptEncoder(nn.Module):
    """Maps input features to interpretable concept scores via sigmoid."""

    def __init__(self, n_input, n_concepts, hidden_dim=256):
        super().__init__()
        self.encoder = nn.Sequential(
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
        """Returns concept activations in [0, 1]."""
        return torch.sigmoid(self.encoder(x))


class PredictionLayer(nn.Module):
    """Multi-task prediction from concept scores using CORAL ordinal heads."""

    def __init__(self, n_concepts):
        super().__init__()
        hidden = max(64, n_concepts)

        # Shared projection from concepts
        self.shared = nn.Sequential(
            nn.Linear(n_concepts, hidden),
            nn.BatchNorm1d(hidden),
            nn.ReLU(),
            nn.Dropout(0.2),
        )

        # Disease binary head (standard classification)
        self.head_disease = nn.Linear(hidden, 2)

        # CORAL ordinal heads
        self.head_nas_3class = CoralLayer(hidden, 3)
        self.head_nas_9class = CoralLayer(hidden, 9)
        self.head_steatosis = CoralLayer(hidden, 4)
        self.head_inflammation = CoralLayer(hidden, 3)
        self.head_ballooning = CoralLayer(hidden, 3)
        self.head_fibrosis = CoralLayer(hidden, 5)

    def forward(self, concepts):
        z = self.shared(concepts)
        return {
            "disease": self.head_disease(z),
            "nas_3class": self.head_nas_3class(z),
            "nas_9class": self.head_nas_9class(z),
            "steatosis": self.head_steatosis(z),
            "inflammation": self.head_inflammation(z),
            "ballooning": self.head_ballooning(z),
            "fibrosis": self.head_fibrosis(z),
        }


class ConceptBottleneckModel(nn.Module):
    """Full concept bottleneck: encoder -> concepts -> predictions."""

    def __init__(self, n_input, n_concepts, hidden_dim=256):
        super().__init__()
        self.concept_encoder = ConceptEncoder(n_input, n_concepts, hidden_dim)
        self.predictor = PredictionLayer(n_concepts)

    def forward(self, x):
        concepts = self.concept_encoder(x)
        task_logits = self.predictor(concepts)
        return concepts, task_logits


# =============================================================================
# Kendall multi-task uncertainty weighting
# =============================================================================
class KendallUncertaintyWeights(nn.Module):
    """Learnable task-specific uncertainty parameters (Kendall et al., 2018).

    Each task has a log-variance parameter log_sigma^2.
    Loss_total = sum_t (1/(2*sigma_t^2)) * loss_t + log(sigma_t)
    """

    def __init__(self, n_tasks, init_vals=None):
        super().__init__()
        if init_vals is not None:
            self.log_vars = nn.Parameter(
                torch.tensor(init_vals, dtype=torch.float32)
            )
        else:
            self.log_vars = nn.Parameter(torch.zeros(n_tasks))

    def forward(self, losses):
        """
        Args:
            losses: list of scalar tensors (one per task)
        Returns:
            weighted total loss (scalar)
        """
        total = torch.tensor(0.0, device=self.log_vars.device)
        for i, loss_i in enumerate(losses):
            precision = torch.exp(-self.log_vars[i])
            total = total + precision * loss_i + self.log_vars[i]
        return total


# =============================================================================
# Dataset
# =============================================================================
class ConceptDataset(Dataset):
    """Features + per-task labels + concept supervision targets."""

    def __init__(self, X, label_dict, concept_targets=None):
        self.X = torch.tensor(X, dtype=torch.float32)
        self.labels = {
            k: torch.tensor(v, dtype=torch.long)
            for k, v in label_dict.items()
        }
        # Concept targets: dict of name -> tensor (supervised) or None
        self.concept_targets = None
        if concept_targets is not None:
            self.concept_targets = torch.tensor(
                np.column_stack(list(concept_targets.values())),
                dtype=torch.float32,
            )
            # Mask: -1 means unsupervised
            self.concept_mask = (self.concept_targets >= 0).float()

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        item = {
            "X": self.X[idx],
            "labels": {k: v[idx] for k, v in self.labels.items()},
        }
        if self.concept_targets is not None:
            item["concept_targets"] = self.concept_targets[idx]
            item["concept_mask"] = self.concept_mask[idx]
        return item


# =============================================================================
# Evaluation
# =============================================================================
def quadratic_weighted_kappa(y_true, y_pred):
    """QWK implementation."""
    classes = np.union1d(y_true, y_pred).astype(int)
    n_classes = int(classes.max()) + 1
    if n_classes < 2:
        return 0.0
    n = len(y_true)
    cm = np.zeros((n_classes, n_classes), dtype=float)
    for t, p in zip(y_true.astype(int), y_pred.astype(int)):
        cm[t, p] += 1
    w = np.zeros((n_classes, n_classes))
    for i in range(n_classes):
        for j in range(n_classes):
            w[i, j] = (i - j) ** 2 / max((n_classes - 1) ** 2, 1)
    hist_t = cm.sum(axis=1)
    hist_p = cm.sum(axis=0)
    expected = np.outer(hist_t, hist_p) / max(n, 1)
    num = (w * cm).sum()
    den = (w * expected).sum()
    return 1.0 - num / max(den, 1e-10)


def eval_ordinal(y_true, probs, task_name, n_classes):
    """Evaluate ordinal task: QWK, MAE, accuracy, adjacent accuracy."""
    valid = y_true >= 0
    n_valid = int(valid.sum())
    if n_valid < 5:
        return {"task": task_name, "n_valid": n_valid,
                "qwk": np.nan, "mae": np.nan, "acc": np.nan, "adj_acc": np.nan}
    y_v = y_true[valid].astype(int)
    y_pred = probs[valid].argmax(axis=1)
    qwk = quadratic_weighted_kappa(y_v, y_pred)
    mae = float(np.mean(np.abs(y_v - y_pred)))
    acc = float(accuracy_score(y_v, y_pred))
    adj_acc = float(np.mean(np.abs(y_v - y_pred) <= 1))
    return {"task": task_name, "n_valid": n_valid,
            "qwk": qwk, "mae": mae, "acc": acc, "adj_acc": adj_acc}


def eval_binary(y_true, probs, task_name):
    """Evaluate binary classification: AUROC, accuracy."""
    valid = y_true >= 0
    n_valid = int(valid.sum())
    if n_valid < 5:
        return {"task": task_name, "n_valid": n_valid,
                "auroc": np.nan, "acc": np.nan}
    y_v = y_true[valid]
    p_v = probs[valid]
    y_pred = p_v.argmax(axis=1)
    acc = float(accuracy_score(y_v, y_pred))
    try:
        auroc = float(roc_auc_score(y_v, p_v[:, 1]))
    except ValueError:
        auroc = np.nan
    return {"task": task_name, "n_valid": n_valid, "auroc": auroc, "acc": acc}


def compute_class_weights(y, n_classes, dev):
    """Inverse-frequency weights for balanced CE; ignores -1 labels."""
    valid = y[y >= 0]
    if len(valid) == 0:
        return torch.ones(n_classes, device=dev)
    counts = np.bincount(valid, minlength=n_classes).astype(np.float32)
    counts = np.maximum(counts, 1.0)
    weights = len(valid) / (n_classes * counts)
    return torch.tensor(weights, dtype=torch.float32, device=dev)


# =============================================================================
# Training
# =============================================================================
def train_one_fold(X_train, lab_train, concept_train,
                   X_val, lab_val, concept_val,
                   n_features, fold_name="fold",
                   phase1_epochs=100, phase2_epochs=200, finetune_epochs=100,
                   lr=5e-4, batch_size=128, patience=30):
    """Train concept bottleneck for one LOCO fold.

    Three training phases:
      Phase 1: Train concept encoder to predict supervised concept targets
      Phase 2: Freeze concept encoder, train prediction heads
      Phase 3: End-to-end fine-tuning at reduced LR
    """
    log.info(f"\n  [{fold_name}] Training concept bottleneck...")
    log.info(f"    Train: {X_train.shape[0]}, Val: {X_val.shape[0]}")

    train_ds = ConceptDataset(X_train, lab_train, concept_train)
    val_ds = ConceptDataset(X_val, lab_val, concept_val)

    train_loader = DataLoader(
        train_ds, batch_size=batch_size, shuffle=True,
        drop_last=False, num_workers=0, pin_memory=True,
    )
    val_loader = DataLoader(
        val_ds, batch_size=batch_size, shuffle=False,
        num_workers=0, pin_memory=True,
    )

    # Disease class weights
    cw_disease = compute_class_weights(lab_train["disease"], 2, device)

    model = ConceptBottleneckModel(
        n_input=n_features,
        n_concepts=N_CONCEPTS,
        hidden_dim=min(512, max(256, n_features // 2)),
    ).to(device)

    # Kendall uncertainty weights (one per task)
    task_names = list(TASK_CFG.keys())
    uncertainty = KendallUncertaintyWeights(len(task_names)).to(device)

    history = []

    # ------ Phase 1: Concept Supervision ------
    log.info(f"  [{fold_name}] Phase 1: Concept supervision ({phase1_epochs} epochs)")
    concept_optimizer = optim.Adam(
        model.concept_encoder.parameters(), lr=lr, weight_decay=1e-4
    )
    concept_scheduler = optim.lr_scheduler.CosineAnnealingLR(
        concept_optimizer, T_max=phase1_epochs, eta_min=1e-6
    )

    best_concept_loss = float("inf")
    best_concept_state = None
    p1_patience = 0

    for epoch in range(phase1_epochs):
        model.concept_encoder.train()
        total_loss = 0.0
        n_batches = 0

        for batch in train_loader:
            X_b = batch["X"].to(device, non_blocking=True)
            ct = batch["concept_targets"].to(device, non_blocking=True)
            cm = batch["concept_mask"].to(device, non_blocking=True)

            concepts = model.concept_encoder(X_b)

            # MSE loss on supervised concepts only
            diff = (concepts - ct) ** 2
            masked_diff = diff * cm
            if cm.sum() > 0:
                loss = masked_diff.sum() / cm.sum()
            else:
                loss = torch.tensor(0.0, device=device, requires_grad=True)

            concept_optimizer.zero_grad()
            loss.backward()
            nn.utils.clip_grad_norm_(model.concept_encoder.parameters(), 5.0)
            concept_optimizer.step()

            total_loss += loss.item()
            n_batches += 1

        concept_scheduler.step()
        avg_loss = total_loss / max(n_batches, 1)

        # Validate concept reconstruction
        model.concept_encoder.eval()
        val_closs = 0.0
        n_vb = 0
        with torch.no_grad():
            for batch in val_loader:
                X_b = batch["X"].to(device, non_blocking=True)
                ct = batch["concept_targets"].to(device, non_blocking=True)
                cm = batch["concept_mask"].to(device, non_blocking=True)
                concepts = model.concept_encoder(X_b)
                diff = (concepts - ct) ** 2
                masked_diff = diff * cm
                if cm.sum() > 0:
                    val_closs += (masked_diff.sum() / cm.sum()).item()
                n_vb += 1

        val_closs /= max(n_vb, 1)
        history.append({
            "fold": fold_name, "phase": 1, "epoch": epoch + 1,
            "train_concept_loss": avg_loss, "val_concept_loss": val_closs,
        })

        if val_closs < best_concept_loss:
            best_concept_loss = val_closs
            best_concept_state = {
                k: v.cpu().clone()
                for k, v in model.concept_encoder.state_dict().items()
            }
            p1_patience = 0
        else:
            p1_patience += 1

        if p1_patience >= patience:
            log.info(f"    Phase 1 early stop at epoch {epoch + 1}")
            break

        if (epoch + 1) % 25 == 0 or epoch == 0:
            log.info(f"    P1 e={epoch+1}: train={avg_loss:.4f} val={val_closs:.4f}")

    # Load best concept encoder
    if best_concept_state is not None:
        model.concept_encoder.load_state_dict(best_concept_state)
    model.concept_encoder.to(device)

    # ------ Phase 2: Prediction Training (frozen concepts) ------
    log.info(f"  [{fold_name}] Phase 2: Prediction training ({phase2_epochs} epochs)")
    for p in model.concept_encoder.parameters():
        p.requires_grad = False

    pred_params = list(model.predictor.parameters()) + list(uncertainty.parameters())
    pred_optimizer = optim.Adam(pred_params, lr=lr, weight_decay=1e-4)
    pred_scheduler = optim.lr_scheduler.CosineAnnealingLR(
        pred_optimizer, T_max=phase2_epochs, eta_min=1e-6
    )

    best_val_loss = float("inf")
    best_pred_state = None
    best_unc_state = None
    p2_patience = 0

    for epoch in range(phase2_epochs):
        model.predictor.train()
        model.concept_encoder.eval()
        train_losses = defaultdict(float)
        n_batches = 0

        for batch in train_loader:
            X_b = batch["X"].to(device, non_blocking=True)
            y_b = {k: v.to(device, non_blocking=True)
                   for k, v in batch["labels"].items()}

            with torch.no_grad():
                concepts = model.concept_encoder(X_b)
            task_logits = model.predictor(concepts)

            # Per-task losses
            task_losses = []

            # Disease: CE
            tl = masked_cross_entropy(
                task_logits["disease"], y_b["disease"], weight=cw_disease
            )
            task_losses.append(tl)
            train_losses["disease"] += tl.item()

            # Ordinal tasks: CORAL loss
            for task in ["nas_3class", "nas_9class", "steatosis",
                         "inflammation", "ballooning", "fibrosis"]:
                n_cls = TASK_CFG[task][0]
                tl = masked_coral_loss(task_logits[task], y_b[task], n_cls)
                task_losses.append(tl)
                train_losses[task] += tl.item()

            # Kendall uncertainty-weighted total
            loss = uncertainty(task_losses)

            pred_optimizer.zero_grad()
            loss.backward()
            nn.utils.clip_grad_norm_(pred_params, 5.0)
            pred_optimizer.step()
            n_batches += 1

        pred_scheduler.step()

        # Validate
        model.eval()
        val_losses = defaultdict(float)
        n_vb = 0
        with torch.no_grad():
            for batch in val_loader:
                X_b = batch["X"].to(device, non_blocking=True)
                y_b = {k: v.to(device, non_blocking=True)
                       for k, v in batch["labels"].items()}
                concepts = model.concept_encoder(X_b)
                task_logits = model.predictor(concepts)

                vl = masked_cross_entropy(
                    task_logits["disease"], y_b["disease"], weight=cw_disease
                )
                val_losses["disease"] += vl.item()
                for task in ["nas_3class", "nas_9class", "steatosis",
                             "inflammation", "ballooning", "fibrosis"]:
                    n_cls = TASK_CFG[task][0]
                    vl = masked_coral_loss(task_logits[task], y_b[task], n_cls)
                    val_losses[task] += vl.item()
                n_vb += 1

        total_val = sum(
            TASK_CFG[t][1] * val_losses[t] / max(n_vb, 1)
            for t in TASK_CFG
        )

        row = {
            "fold": fold_name, "phase": 2, "epoch": epoch + 1,
            "total_val_loss": total_val,
        }
        for t in TASK_CFG:
            row[f"train_{t}"] = train_losses[t] / max(n_batches, 1)
            row[f"val_{t}"] = val_losses[t] / max(n_vb, 1)
        # Record uncertainty weights
        for i, t in enumerate(task_names):
            row[f"sigma_{t}"] = float(torch.exp(0.5 * uncertainty.log_vars[i]).item())
        history.append(row)

        if total_val < best_val_loss:
            best_val_loss = total_val
            best_pred_state = {
                k: v.cpu().clone()
                for k, v in model.predictor.state_dict().items()
            }
            best_unc_state = {
                k: v.cpu().clone()
                for k, v in uncertainty.state_dict().items()
            }
            p2_patience = 0
        else:
            p2_patience += 1

        if p2_patience >= patience:
            log.info(f"    Phase 2 early stop at epoch {epoch + 1}")
            break

        if (epoch + 1) % 25 == 0 or epoch == 0:
            parts = [f"e={epoch+1:3d}"]
            for t in list(TASK_CFG.keys())[:4]:
                parts.append(f"{t[:6]}={train_losses[t]/max(n_batches,1):.3f}")
            parts.append(f"val={total_val:.4f}")
            log.info(f"    P2 {' | '.join(parts)}")

    # Load best prediction state
    if best_pred_state is not None:
        model.predictor.load_state_dict(best_pred_state)
    if best_unc_state is not None:
        uncertainty.load_state_dict(best_unc_state)
    model.predictor.to(device)
    uncertainty.to(device)

    # ------ Phase 3: End-to-End Fine-tuning ------
    log.info(f"  [{fold_name}] Phase 3: End-to-end fine-tuning ({finetune_epochs} epochs)")
    for p in model.concept_encoder.parameters():
        p.requires_grad = True

    all_params = list(model.parameters()) + list(uncertainty.parameters())
    ft_optimizer = optim.Adam(all_params, lr=lr * 0.1, weight_decay=1e-4)
    ft_scheduler = optim.lr_scheduler.CosineAnnealingLR(
        ft_optimizer, T_max=finetune_epochs, eta_min=1e-7
    )

    best_ft_loss = float("inf")
    best_ft_state = None
    p3_patience = 0

    for epoch in range(finetune_epochs):
        model.train()
        train_losses = defaultdict(float)
        n_batches = 0

        for batch in train_loader:
            X_b = batch["X"].to(device, non_blocking=True)
            y_b = {k: v.to(device, non_blocking=True)
                   for k, v in batch["labels"].items()}
            ct = batch["concept_targets"].to(device, non_blocking=True)
            cm = batch["concept_mask"].to(device, non_blocking=True)

            concepts, task_logits = model(X_b)

            # Task losses
            task_losses = []
            tl = masked_cross_entropy(
                task_logits["disease"], y_b["disease"], weight=cw_disease
            )
            task_losses.append(tl)
            train_losses["disease"] += tl.item()

            for task in ["nas_3class", "nas_9class", "steatosis",
                         "inflammation", "ballooning", "fibrosis"]:
                n_cls = TASK_CFG[task][0]
                tl = masked_coral_loss(task_logits[task], y_b[task], n_cls)
                task_losses.append(tl)
                train_losses[task] += tl.item()

            # Uncertainty-weighted task loss
            loss = uncertainty(task_losses)

            # Concept regularization (keep concepts aligned during fine-tuning)
            diff = (concepts - ct) ** 2
            masked_diff = diff * cm
            if cm.sum() > 0:
                concept_reg = masked_diff.sum() / cm.sum() * 0.1
                loss = loss + concept_reg
                train_losses["concept_reg"] += concept_reg.item()

            ft_optimizer.zero_grad()
            loss.backward()
            nn.utils.clip_grad_norm_(all_params, 5.0)
            ft_optimizer.step()
            n_batches += 1

        ft_scheduler.step()

        # Validate
        model.eval()
        val_losses = defaultdict(float)
        n_vb = 0
        with torch.no_grad():
            for batch in val_loader:
                X_b = batch["X"].to(device, non_blocking=True)
                y_b = {k: v.to(device, non_blocking=True)
                       for k, v in batch["labels"].items()}
                concepts, task_logits = model(X_b)

                vl = masked_cross_entropy(
                    task_logits["disease"], y_b["disease"], weight=cw_disease
                )
                val_losses["disease"] += vl.item()
                for task in ["nas_3class", "nas_9class", "steatosis",
                             "inflammation", "ballooning", "fibrosis"]:
                    n_cls = TASK_CFG[task][0]
                    vl = masked_coral_loss(task_logits[task], y_b[task], n_cls)
                    val_losses[task] += vl.item()
                n_vb += 1

        total_val = sum(
            TASK_CFG[t][1] * val_losses[t] / max(n_vb, 1)
            for t in TASK_CFG
        )

        row = {
            "fold": fold_name, "phase": 3, "epoch": epoch + 1,
            "total_val_loss": total_val,
        }
        for t in TASK_CFG:
            row[f"train_{t}"] = train_losses[t] / max(n_batches, 1)
            row[f"val_{t}"] = val_losses[t] / max(n_vb, 1)
        history.append(row)

        if total_val < best_ft_loss:
            best_ft_loss = total_val
            best_ft_state = {
                k: v.cpu().clone()
                for k, v in model.state_dict().items()
            }
            p3_patience = 0
        else:
            p3_patience += 1

        if p3_patience >= patience:
            log.info(f"    Phase 3 early stop at epoch {epoch + 1}")
            break

        if (epoch + 1) % 25 == 0 or epoch == 0:
            parts = [f"e={epoch+1:3d}"]
            for t in list(TASK_CFG.keys())[:4]:
                parts.append(f"{t[:6]}={train_losses[t]/max(n_batches,1):.3f}")
            parts.append(f"val={total_val:.4f}")
            log.info(f"    P3 {' | '.join(parts)}")

    # Load best fine-tuned state
    if best_ft_state is not None:
        model.load_state_dict(best_ft_state)
    model.to(device)

    # ------ Generate predictions + concept activations on val ------
    model.eval()
    all_logits = {t: [] for t in TASK_CFG}
    all_concepts = []

    with torch.no_grad():
        for batch in val_loader:
            X_b = batch["X"].to(device, non_blocking=True)
            concepts, task_logits = model(X_b)
            all_concepts.append(concepts.cpu().numpy())
            for t in TASK_CFG:
                all_logits[t].append(task_logits[t].cpu())

    # Convert logits to probabilities
    predictions = {}
    for t in TASK_CFG:
        cat = torch.cat(all_logits[t], dim=0)
        n_cls, _, is_ord = TASK_CFG[t]
        if is_ord:
            probs = coral_ordinal_probs(cat, n_cls).numpy()
        else:
            probs = torch.softmax(cat, dim=1).numpy()
        predictions[t] = probs

    val_concepts = np.vstack(all_concepts)

    return predictions, val_concepts, pd.DataFrame(history), model


# =============================================================================
# Main training loop with LOCO cross-validation
# =============================================================================
def main():
    t0 = time.time()

    if not os.path.exists(H5_PATH):
        sys.exit(f"ERROR: {H5_PATH} not found. Run Script 62 first.")
    if not os.path.exists(META_PATH):
        sys.exit(f"ERROR: {META_PATH} not found. Run Script 89 first.")

    if device.type == "cuda":
        log.info(f"  GPU: {torch.cuda.get_device_name(0)}")
        log.info(f"  VRAM: {torch.cuda.get_device_properties(0).total_mem / 1e9:.1f} GB")

    # Load all features and labels
    X, labels, sample_ids, datasets, feature_groups, concept_targets = load_features()
    n_samples, n_features = X.shape

    # ===== NAS LOCO (5-fold) =====
    log.info("\n" + "=" * 60)
    log.info("NAS LOCO Cross-Validation (5 folds)")
    log.info("=" * 60)

    all_results = []
    all_concepts_nas = np.zeros((n_samples, N_CONCEPTS))
    concept_assigned = np.zeros(n_samples, dtype=bool)
    all_history = []

    for fold_ds in NAS_DATASETS:
        fold_name = f"nas_{fold_ds}"
        val_mask = datasets == fold_ds
        train_mask = ~val_mask
        n_val = val_mask.sum()

        if n_val == 0:
            log.warning(f"  [{fold_name}] No validation samples, skipping")
            continue

        # Split data
        X_train, X_val = X[train_mask], X[val_mask]
        lab_train = {t: v[train_mask] for t, v in labels.items()}
        lab_val = {t: v[val_mask] for t, v in labels.items()}

        # Normalize concept targets using train-fold statistics only
        normed_targets = normalize_concept_targets(concept_targets, train_mask)
        ct_train = {k: v[train_mask] for k, v in normed_targets.items()}
        ct_val = {k: v[val_mask] for k, v in normed_targets.items()}

        # Train
        predictions, val_concepts, history, model = train_one_fold(
            X_train, lab_train, ct_train,
            X_val, lab_val, ct_val,
            n_features=n_features,
            fold_name=fold_name,
        )
        all_history.append(history)

        # Store concept activations
        val_indices = np.where(val_mask)[0]
        all_concepts_nas[val_indices] = val_concepts
        concept_assigned[val_indices] = True

        # Evaluate
        for task in TASK_CFG:
            n_cls, _, is_ord = TASK_CFG[task]
            if is_ord:
                metrics = eval_ordinal(
                    lab_val[task], predictions[task], task, n_cls
                )
            else:
                metrics = eval_binary(
                    lab_val[task], predictions[task], task
                )
            metrics["fold"] = fold_name
            metrics["fold_type"] = "nas_loco"
            all_results.append(metrics)

            if is_ord and not np.isnan(metrics.get("qwk", np.nan)):
                log.info(
                    f"  [{fold_name}] {task}: QWK={metrics['qwk']:.3f} "
                    f"MAE={metrics['mae']:.2f} Acc={metrics['acc']:.3f}"
                )
            elif not is_ord and not np.isnan(metrics.get("auroc", np.nan)):
                log.info(
                    f"  [{fold_name}] {task}: AUROC={metrics['auroc']:.3f} "
                    f"Acc={metrics['acc']:.3f}"
                )

    # ===== Fibrosis LOCO (6-fold) =====
    log.info("\n" + "=" * 60)
    log.info("Fibrosis LOCO Cross-Validation (6 folds)")
    log.info("=" * 60)

    for fold_ds in FIB_DATASETS:
        fold_name = f"fib_{fold_ds}"
        val_mask = datasets == fold_ds
        train_mask = ~val_mask
        n_val = val_mask.sum()

        if n_val == 0:
            log.warning(f"  [{fold_name}] No validation samples, skipping")
            continue

        X_train, X_val = X[train_mask], X[val_mask]
        lab_train = {t: v[train_mask] for t, v in labels.items()}
        lab_val = {t: v[val_mask] for t, v in labels.items()}

        # Normalize concept targets using train-fold statistics only
        normed_targets = normalize_concept_targets(concept_targets, train_mask)
        ct_train = {k: v[train_mask] for k, v in normed_targets.items()}
        ct_val = {k: v[val_mask] for k, v in normed_targets.items()}

        predictions, val_concepts, history, model = train_one_fold(
            X_train, lab_train, ct_train,
            X_val, lab_val, ct_val,
            n_features=n_features,
            fold_name=fold_name,
        )
        all_history.append(history)

        # Store concept activations for unassigned samples
        val_indices = np.where(val_mask)[0]
        unassigned = ~concept_assigned[val_indices]
        if unassigned.any():
            all_concepts_nas[val_indices[unassigned]] = val_concepts[unassigned]
            concept_assigned[val_indices[unassigned]] = True

        for task in TASK_CFG:
            n_cls, _, is_ord = TASK_CFG[task]
            if is_ord:
                metrics = eval_ordinal(
                    lab_val[task], predictions[task], task, n_cls
                )
            else:
                metrics = eval_binary(
                    lab_val[task], predictions[task], task
                )
            metrics["fold"] = fold_name
            metrics["fold_type"] = "fib_loco"
            all_results.append(metrics)

            if is_ord and not np.isnan(metrics.get("qwk", np.nan)):
                log.info(
                    f"  [{fold_name}] {task}: QWK={metrics['qwk']:.3f} "
                    f"MAE={metrics['mae']:.2f} Acc={metrics['acc']:.3f}"
                )

    # ===== Train full model for concept activations + save =====
    log.info("\n" + "=" * 60)
    log.info("Training full model (all data) for concept export")
    log.info("=" * 60)

    full_ds = ConceptDataset(X, labels, concept_targets)
    full_loader = DataLoader(
        full_ds, batch_size=128, shuffle=True,
        num_workers=0, pin_memory=True,
    )

    full_model = ConceptBottleneckModel(
        n_input=n_features,
        n_concepts=N_CONCEPTS,
        hidden_dim=min(512, max(256, n_features // 2)),
    ).to(device)
    full_unc = KendallUncertaintyWeights(len(TASK_CFG)).to(device)
    cw_disease = compute_class_weights(labels["disease"], 2, device)

    # Simplified training for full model (single phase, 200 epochs)
    all_params = list(full_model.parameters()) + list(full_unc.parameters())
    optimizer = optim.Adam(all_params, lr=5e-4, weight_decay=1e-4)
    scheduler = optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=200, eta_min=1e-6
    )

    best_loss = float("inf")
    best_state = None
    no_improve = 0

    for epoch in range(200):
        full_model.train()
        total_loss = 0.0
        n_b = 0

        for batch in full_loader:
            X_b = batch["X"].to(device, non_blocking=True)
            y_b = {k: v.to(device, non_blocking=True)
                   for k, v in batch["labels"].items()}
            ct = batch["concept_targets"].to(device, non_blocking=True)
            cm = batch["concept_mask"].to(device, non_blocking=True)

            concepts, task_logits = full_model(X_b)

            task_losses = []
            tl = masked_cross_entropy(
                task_logits["disease"], y_b["disease"], weight=cw_disease
            )
            task_losses.append(tl)
            for task in ["nas_3class", "nas_9class", "steatosis",
                         "inflammation", "ballooning", "fibrosis"]:
                n_cls = TASK_CFG[task][0]
                tl = masked_coral_loss(task_logits[task], y_b[task], n_cls)
                task_losses.append(tl)

            loss = full_unc(task_losses)

            # Concept regularization
            diff = (concepts - ct) ** 2
            masked_diff = diff * cm
            if cm.sum() > 0:
                loss = loss + 0.1 * masked_diff.sum() / cm.sum()

            optimizer.zero_grad()
            loss.backward()
            nn.utils.clip_grad_norm_(all_params, 5.0)
            optimizer.step()
            total_loss += loss.item()
            n_b += 1

        scheduler.step()
        avg_loss = total_loss / max(n_b, 1)

        if avg_loss < best_loss:
            best_loss = avg_loss
            best_state = {
                k: v.cpu().clone()
                for k, v in full_model.state_dict().items()
            }
            no_improve = 0
        else:
            no_improve += 1

        if no_improve >= 50:
            log.info(f"  Full model early stop at epoch {epoch + 1}")
            break

        if (epoch + 1) % 50 == 0 or epoch == 0:
            log.info(f"  Full model e={epoch+1}: loss={avg_loss:.4f}")

    # Load best full model
    if best_state is not None:
        full_model.load_state_dict(best_state)
    full_model.to(device)

    # Extract concept activations for all samples
    full_model.eval()
    all_concepts_full = []
    with torch.no_grad():
        full_eval_loader = DataLoader(
            full_ds, batch_size=128, shuffle=False,
            num_workers=0, pin_memory=True,
        )
        for batch in full_eval_loader:
            X_b = batch["X"].to(device, non_blocking=True)
            concepts, _ = full_model(X_b)
            all_concepts_full.append(concepts.cpu().numpy())

    concepts_full = np.vstack(all_concepts_full)

    # Save model
    torch.save(full_model.state_dict(), OUT_MODEL)
    log.info(f"  Model saved: {OUT_MODEL}")

    # ===== Save results =====
    log.info("\n" + "=" * 60)
    log.info("Saving results")
    log.info("=" * 60)

    # 1. Per-fold per-task results
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(OUT_RESULTS, index=False)
    log.info(f"  Results: {OUT_RESULTS} ({len(results_df)} rows)")

    # 2. Aggregated summary
    summary_rows = []
    for fold_type in ["nas_loco", "fib_loco"]:
        sub = results_df[results_df["fold_type"] == fold_type]
        for task in TASK_CFG:
            task_sub = sub[sub["task"] == task]
            if len(task_sub) == 0:
                continue
            row = {
                "fold_type": fold_type,
                "task": task,
                "n_folds": len(task_sub),
            }
            for metric in ["qwk", "mae", "acc", "adj_acc", "auroc"]:
                if metric in task_sub.columns:
                    vals = task_sub[metric].dropna()
                    if len(vals) > 0:
                        row[f"{metric}_mean"] = vals.mean()
                        row[f"{metric}_std"] = vals.std()
            summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT_SUMMARY, index=False)
    log.info(f"  Summary: {OUT_SUMMARY}")

    # 3. Concept activations (from full model)
    concept_names = list(concept_targets.keys())
    concept_df = pd.DataFrame(concepts_full, columns=concept_names)
    concept_df.insert(0, "sample_id", sample_ids)
    concept_df.to_csv(OUT_CONCEPTS, index=False)
    log.info(f"  Concepts: {OUT_CONCEPTS} ({concept_df.shape})")

    # 4. Training history
    history_df = pd.concat(all_history, ignore_index=True)
    history_path = os.path.join(OUTDIR, "concept_bottleneck_training_curves.csv")
    history_df.to_csv(history_path, index=False)

    # Print summary
    log.info("\n" + "=" * 60)
    log.info("SUMMARY")
    log.info("=" * 60)
    log.info(f"  Features: {n_features} input -> {N_CONCEPTS} concepts -> 7 tasks")
    for _, row in summary_df.iterrows():
        ft = row["fold_type"]
        task = row["task"]
        if "qwk_mean" in row and pd.notna(row.get("qwk_mean")):
            log.info(
                f"  {ft} {task}: QWK={row['qwk_mean']:.3f}+/-{row.get('qwk_std',0):.3f}"
            )
        elif "auroc_mean" in row and pd.notna(row.get("auroc_mean")):
            log.info(
                f"  {ft} {task}: AUROC={row['auroc_mean']:.3f}+/-{row.get('auroc_std',0):.3f}"
            )

    elapsed = time.time() - t0
    log.info(f"\nDone in {elapsed/60:.1f} minutes.")
    log.info(f"Outputs:")
    log.info(f"  {OUT_RESULTS}")
    log.info(f"  {OUT_SUMMARY}")
    log.info(f"  {OUT_CONCEPTS}")
    log.info(f"  {OUT_MODEL}")


if __name__ == "__main__":
    main()
