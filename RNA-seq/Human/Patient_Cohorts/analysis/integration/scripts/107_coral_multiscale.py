#!/usr/bin/env python3
"""
107_coral_multiscale.py
CORAL ordinal neural network on multi-scale features with proper LOCO-CV.

Architecture:
  CoralMLP: shared trunk (Linear 256 -> BN -> ReLU -> Drop 0.3 ->
            Linear 64 -> BN -> ReLU -> Drop 0.2) -> CoralLayer.

Features per sample (~149 total):
  - NAS-VAE 64-dim LOCO embeddings (from nas_embeddings_loco/)
  - ssGSEA pathway scores: 1791 pathways PCA'd to 50 dims (pathway_scores_ssgsea.rds)
  - TF activities: 295 TFs PCA'd to 30 dims (tf_activity_features.rds)
  - Deconv proportions: 4 BayesPrism columns (unified_bayesprism_proportions.csv)
  - Sex: 1 binary feature
  Total: 64 + 50 + 30 + 4 + 1 = 149

Training: 5-fold NAS LOCO + 6-fold fibrosis LOCO. For each fold:
  1. Load LOCO-specific NAS embeddings (model trained WITHOUT held-out cohort)
  2. Assemble multi-scale feature matrix
  3. PCA fitted on training fold only (no leakage)
  4. Train CoralMLP with coral_loss for 200 epochs, early stopping patience=20
  5. Evaluate on held-out fold: QWK, AUROC, MAE

Targets:
  - NAS 9-class (0-8)
  - NAS 4-group
  - NAS>=5 binary
  - Fibrosis F0-F4
  - F>=3 binary

Input (all from results/staging_classifier/):
  - nas_embeddings_loco/nas_embeddings_loco_{fold}.csv
  - pathway_scores_ssgsea.rds   (1791 x 1444, pathways x samples — TRANSPOSE)
  - tf_activity_features.rds    (1444 x 295, samples x TFs)
  - unified_bayesprism_proportions.csv (deconv from Script 22)
  - modeling_metadata.csv        (labels + LOCO fold assignments)

Output (all to results/staging_classifier/):
  - coral_multiscale_results.csv   (per-fold per-target metrics)
  - coral_multiscale_summary.csv   (aggregated mean +/- std)

SLURM: cpu partition, 8 CPUs, 64G RAM, 48h
Env:   micromamba activate rapids_singlecell

Usage:
  sbatch --job-name=stg107_coral \
         --partition=cpu --cpus-per-task=8 --mem=64G --time=48:00:00 \
         --output=logs/107_coral_%j.out \
         --error=logs/107_coral_%j.err \
         --wrap="bash -c 'eval \"\\$(micromamba shell hook --shell bash)\" && \\
                 micromamba activate rapids_singlecell && \\
                 cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts && \\
                 python 107_coral_multiscale.py'"
"""

# Must be set before ANY other imports to prevent CUDA crash on cpu nodes
import os
os.environ["RAPIDS_NO_INITIALIZE"] = "1"
os.environ["CUDF_PANDAS"] = "0"

import sys
import time
import subprocess
import warnings
import logging
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader

from coral_pytorch.layers import CoralLayer
from coral_pytorch.losses import coral_loss as coral_loss_fn
from coral_pytorch.dataset import levels_from_labelbatch, proba_to_label

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    roc_auc_score, mean_absolute_error, confusion_matrix
)

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
META_PATH = os.path.join(OUTDIR, "modeling_metadata.csv")
DECONV_PATH = os.path.join(
    INT, "results/deconvolution/bayesprism/unified_bayesprism_proportions.csv"
)
LOCO_EMB_DIR = os.path.join(OUTDIR, "nas_embeddings_loco")

# RDS inputs (from upstream scripts 103/105)
SSGSEA_RDS = os.path.join(OUTDIR, "pathway_scores_ssgsea.rds")
TF_RDS = os.path.join(OUTDIR, "tf_activity_features.rds")

# Outputs
OUT_RESULTS = os.path.join(OUTDIR, "coral_multiscale_results.csv")
OUT_SUMMARY = os.path.join(OUTDIR, "coral_multiscale_summary.csv")

# Rscript binary
RSCRIPT = "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"

MISSING = -1

# Device (CPU expected for this script)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(os.path.join(LOGDIR, "107_coral_multiscale.log")),
    ],
)
log = logging.getLogger(__name__)

print("=" * 70)
print("107: CORAL Ordinal Multi-Scale Model for MASLD Staging")
print("=" * 70)
print(f"Device : {device}")
print(f"Output : {OUTDIR}")
print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print()

# ---------------------------------------------------------------------------
# Hyperparameters
# ---------------------------------------------------------------------------
N_PCA_PATHWAY = 50       # PCA dims for 1791 ssGSEA pathways
N_PCA_TF = 30            # PCA dims for 295 TF activities
MAX_EPOCHS = 200
PATIENCE = 20
LR = 1e-3
WEIGHT_DECAY = 1e-4
BATCH_SIZE = 64

# LOCO fold definitions
NAS_DATASETS = ["GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066"]
FIB_DATASETS = [
    "GSE130970", "GSE135251", "GSE162694", "GSE174478",
    "GSE193066", "GSE240729",
]


# =============================================================================
# RDS conversion utility
# =============================================================================
def rds_to_csv(rds_path, csv_path, transpose=False, force=False):
    """Convert an RDS file (matrix/data.frame) to CSV via Rscript."""
    if os.path.exists(csv_path) and not force:
        log.info(f"  CSV cache exists: {csv_path}")
        return csv_path

    if not os.path.exists(rds_path):
        log.warning(f"  RDS not found: {rds_path}")
        return None

    transpose_line = "obj <- t(obj)" if transpose else ""

    r_code = f"""
obj <- readRDS("{rds_path}")
if (is.matrix(obj)) {{
    {transpose_line}
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
    if os.path.exists(r_script):
        os.remove(r_script)
    return csv_path


# =============================================================================
# QWK helper
# =============================================================================
def quadratic_weighted_kappa(y_true, y_pred, n_classes):
    """Compute QWK from integer labels."""
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


# =============================================================================
# Model definition
# =============================================================================
class CoralMLP(nn.Module):
    """MLP with CORAL ordinal output layer."""

    def __init__(self, n_features, n_classes):
        super().__init__()
        self.shared = nn.Sequential(
            nn.Linear(n_features, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(256, 64),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.Dropout(0.2),
        )
        self.coral = CoralLayer(64, n_classes)

    def forward(self, x):
        h = self.shared(x)
        return self.coral(h)


class OrdinalDataset(Dataset):
    """Simple dataset for features + ordinal labels."""

    def __init__(self, X, y):
        self.X = torch.tensor(X, dtype=torch.float32)
        self.y = torch.tensor(y, dtype=torch.long)

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]


# =============================================================================
# Training loop
# =============================================================================
def train_coral_model(X_train, y_train, X_val, y_val, n_classes, target_name,
                      fold_name):
    """
    Train CoralMLP with early stopping.

    Returns dict with predictions and metrics on the validation set.
    """
    n_features = X_train.shape[1]

    # Build data loaders
    train_ds = OrdinalDataset(X_train, y_train)
    val_ds = OrdinalDataset(X_val, y_val)
    train_dl = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True,
                          drop_last=False)
    val_dl = DataLoader(val_ds, batch_size=len(val_ds), shuffle=False)

    model = CoralMLP(n_features, n_classes).to(device)
    optimizer = optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=7, min_lr=1e-6
    )

    best_val_loss = float("inf")
    best_state = None
    epochs_no_improve = 0

    for epoch in range(MAX_EPOCHS):
        # --- Train ---
        model.train()
        train_losses = []
        for xb, yb in train_dl:
            xb, yb = xb.to(device), yb.to(device)
            logits = model(xb)
            levels = levels_from_labelbatch(yb, num_classes=n_classes).to(device)
            loss = coral_loss_fn(logits, levels)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            train_losses.append(loss.item())

        # --- Validate ---
        model.eval()
        with torch.no_grad():
            for xb, yb in val_dl:
                xb, yb = xb.to(device), yb.to(device)
                logits = model(xb)
                levels = levels_from_labelbatch(yb, num_classes=n_classes).to(device)
                val_loss = coral_loss_fn(logits, levels).item()

        scheduler.step(val_loss)

        # Early stopping
        if val_loss < best_val_loss - 1e-5:
            best_val_loss = val_loss
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
            epochs_no_improve = 0
        else:
            epochs_no_improve += 1

        if epochs_no_improve >= PATIENCE:
            log.info(f"    [{target_name}/{fold_name}] Early stop at epoch "
                     f"{epoch + 1} (best loss={best_val_loss:.4f})")
            break

    if epoch == MAX_EPOCHS - 1:
        log.info(f"    [{target_name}/{fold_name}] Ran full {MAX_EPOCHS} epochs "
                 f"(best loss={best_val_loss:.4f})")

    # Reload best model
    if best_state is not None:
        model.load_state_dict(best_state)

    # --- Predict on validation set ---
    model.eval()
    with torch.no_grad():
        for xb, yb in val_dl:
            xb = xb.to(device)
            logits = model(xb)
            # CORAL: cumulative probabilities -> predicted labels
            cum_probs = torch.sigmoid(logits)  # (N, n_classes - 1)
            preds = proba_to_label(cum_probs).cpu().numpy()
            # Class probabilities for AUROC: P(Y=k) via differencing
            cum_probs_cpu = cum_probs.cpu()
            ones = torch.ones(cum_probs_cpu.shape[0], 1)
            zeros = torch.zeros(cum_probs_cpu.shape[0], 1)
            extended = torch.cat([ones, cum_probs_cpu, zeros], dim=1)
            class_probs = extended[:, :-1] - extended[:, 1:]
            class_probs = torch.clamp(class_probs, min=1e-7)  # numerical safety
            # Renormalise to sum to 1
            class_probs = class_probs / class_probs.sum(dim=1, keepdim=True)
            probs_np = class_probs.numpy()

    y_val_np = y_val.copy()

    # --- Metrics ---
    qwk = quadratic_weighted_kappa(y_val_np, preds, n_classes)
    mae = mean_absolute_error(y_val_np, preds)

    # AUROC: macro OVR for multi-class, standard for binary
    auroc = np.nan
    if n_classes == 2:
        try:
            auroc = roc_auc_score(y_val_np, probs_np[:, 1])
        except ValueError:
            pass
    else:
        try:
            auroc = roc_auc_score(
                y_val_np, probs_np, multi_class="ovr", average="macro",
                labels=list(range(n_classes)),
            )
        except ValueError:
            pass

    return {
        "qwk": qwk,
        "auroc": auroc,
        "mae": mae,
        "preds": preds,
        "probs": probs_np,
        "y_true": y_val_np,
    }


# =============================================================================
# Load global feature matrices (ssGSEA, TF, deconv)
# =============================================================================
def load_global_features(sample_ids):
    """
    Load ssGSEA pathway scores, TF activities, and deconv proportions.
    Returns raw DataFrames aligned to sample_ids (before PCA).
    """
    # --- ssGSEA pathway scores (1791 pathways x 1444 samples -> transpose) ---
    ssgsea_csv = os.path.join(OUTDIR, "pathway_scores_ssgsea_transposed.csv")
    ssgsea_df = None
    csv_path = rds_to_csv(SSGSEA_RDS, ssgsea_csv, transpose=True)
    if csv_path and os.path.exists(csv_path):
        log.info("  Loading ssGSEA pathway scores (transposed)...")
        raw = pd.read_csv(csv_path)
        if "sample_id" in raw.columns:
            raw = raw.set_index("sample_id")
        ssgsea_df = raw.reindex(sample_ids).fillna(0.0)
        log.info(f"    ssGSEA: {ssgsea_df.shape[1]} pathways x "
                 f"{ssgsea_df.shape[0]} samples")
    else:
        log.warning("  ssGSEA scores not available — feature block will be empty")
        ssgsea_df = pd.DataFrame(index=sample_ids)

    # --- TF activities (1444 samples x 295 TFs, no transpose needed) ---
    tf_csv = os.path.join(OUTDIR, "tf_activity_features_raw.csv")
    tf_df = None
    csv_path = rds_to_csv(TF_RDS, tf_csv, transpose=False)
    if csv_path and os.path.exists(csv_path):
        log.info("  Loading TF activity features...")
        raw = pd.read_csv(csv_path)
        if "sample_id" in raw.columns:
            raw = raw.set_index("sample_id")
        tf_df = raw.reindex(sample_ids).fillna(0.0)
        log.info(f"    TF: {tf_df.shape[1]} TFs x {tf_df.shape[0]} samples")
    else:
        log.warning("  TF activity features not available — feature block will be "
                     "empty")
        tf_df = pd.DataFrame(index=sample_ids)

    # --- Deconv proportions ---
    deconv_df = None
    if os.path.exists(DECONV_PATH):
        log.info("  Loading deconv proportions...")
        raw = pd.read_csv(DECONV_PATH)
        bp_cols = [c for c in raw.columns if c.startswith("bp_")]
        if bp_cols:
            deconv_df = raw.set_index("sample_id")[bp_cols].reindex(
                sample_ids
            ).fillna(0.0)
            log.info(f"    Deconv: {len(bp_cols)} cell types ({bp_cols})")
    if deconv_df is None:
        log.warning("  Deconv proportions not available — feature block will be "
                     "empty")
        deconv_df = pd.DataFrame(index=sample_ids)

    return ssgsea_df, tf_df, deconv_df


# =============================================================================
# Assemble features for one LOCO fold
# =============================================================================
def assemble_fold_features(
    fold_name,
    train_ids,
    test_ids,
    loco_emb_df,
    ssgsea_df,
    tf_df,
    deconv_df,
    meta_df,
):
    """
    Assemble multi-scale features for a single LOCO fold.
    PCA is fit on training data only (no leakage).

    Returns (X_train, X_test) as numpy arrays.
    """
    # --- 1. VAE embeddings from LOCO-specific file ---
    emb_cols = [c for c in loco_emb_df.columns
                if c not in ("sample_id", "is_heldout")]

    X_emb_train = loco_emb_df.loc[
        loco_emb_df.index.isin(train_ids), emb_cols
    ].reindex(train_ids).fillna(0.0).values
    X_emb_test = loco_emb_df.loc[
        loco_emb_df.index.isin(test_ids), emb_cols
    ].reindex(test_ids).fillna(0.0).values

    log.info(f"    VAE embeddings: {X_emb_train.shape[1]} dims")

    # --- 2. ssGSEA pathway scores -> PCA to N_PCA_PATHWAY dims ---
    if ssgsea_df.shape[1] > 0:
        ssgsea_train_raw = ssgsea_df.reindex(train_ids).fillna(0.0).values
        ssgsea_test_raw = ssgsea_df.reindex(test_ids).fillna(0.0).values

        n_pca_pw = min(N_PCA_PATHWAY, ssgsea_train_raw.shape[1],
                       ssgsea_train_raw.shape[0])
        pw_scaler = StandardScaler()
        ssgsea_train_sc = pw_scaler.fit_transform(ssgsea_train_raw)
        ssgsea_test_sc = pw_scaler.transform(ssgsea_test_raw)

        pw_pca = PCA(n_components=n_pca_pw, random_state=SEED)
        X_pw_train = pw_pca.fit_transform(ssgsea_train_sc)
        X_pw_test = pw_pca.transform(ssgsea_test_sc)
        log.info(f"    ssGSEA PCA: {ssgsea_df.shape[1]} -> {n_pca_pw} dims "
                 f"({pw_pca.explained_variance_ratio_.sum():.1%} var)")
    else:
        X_pw_train = np.zeros((len(train_ids), 0))
        X_pw_test = np.zeros((len(test_ids), 0))

    # --- 3. TF activities -> PCA to N_PCA_TF dims ---
    if tf_df.shape[1] > 0:
        tf_train_raw = tf_df.reindex(train_ids).fillna(0.0).values
        tf_test_raw = tf_df.reindex(test_ids).fillna(0.0).values

        n_pca_tf = min(N_PCA_TF, tf_train_raw.shape[1], tf_train_raw.shape[0])
        tf_scaler = StandardScaler()
        tf_train_sc = tf_scaler.fit_transform(tf_train_raw)
        tf_test_sc = tf_scaler.transform(tf_test_raw)

        tf_pca = PCA(n_components=n_pca_tf, random_state=SEED)
        X_tf_train = tf_pca.fit_transform(tf_train_sc)
        X_tf_test = tf_pca.transform(tf_test_sc)
        log.info(f"    TF PCA: {tf_df.shape[1]} -> {n_pca_tf} dims "
                 f"({tf_pca.explained_variance_ratio_.sum():.1%} var)")
    else:
        X_tf_train = np.zeros((len(train_ids), 0))
        X_tf_test = np.zeros((len(test_ids), 0))

    # --- 4. Deconv proportions (low-dim, no PCA needed) ---
    if deconv_df.shape[1] > 0:
        X_dc_train = deconv_df.reindex(train_ids).fillna(0.0).values
        X_dc_test = deconv_df.reindex(test_ids).fillna(0.0).values
        log.info(f"    Deconv: {deconv_df.shape[1]} dims")
    else:
        X_dc_train = np.zeros((len(train_ids), 0))
        X_dc_test = np.zeros((len(test_ids), 0))

    # --- 5. Sex (binary: M=0, F=1) ---
    sex_map = meta_df.set_index("sample_id")["sex"]
    sex_train = sex_map.reindex(train_ids).map({"M": 0, "F": 1}).fillna(0).values.reshape(-1, 1)
    sex_test = sex_map.reindex(test_ids).map({"M": 0, "F": 1}).fillna(0).values.reshape(-1, 1)

    # --- Concatenate ---
    X_train = np.hstack([
        X_emb_train, X_pw_train, X_tf_train, X_dc_train, sex_train
    ]).astype(np.float32)
    X_test = np.hstack([
        X_emb_test, X_pw_test, X_tf_test, X_dc_test, sex_test
    ]).astype(np.float32)

    log.info(f"    Total features: {X_train.shape[1]} "
             f"(train={X_train.shape[0]}, test={X_test.shape[0]})")

    # --- Final standardisation on concatenated features ---
    final_scaler = StandardScaler()
    X_train = final_scaler.fit_transform(X_train)
    X_test = final_scaler.transform(X_test)

    return X_train, X_test


# =============================================================================
# Main
# =============================================================================
def main():
    t0 = time.time()

    # ------------------------------------------------------------------
    # 1. Load metadata
    # ------------------------------------------------------------------
    log.info("Loading metadata...")
    meta = pd.read_csv(META_PATH)
    meta["sample_id"] = meta["sample_id"].astype(str).str.strip()
    all_sample_ids = meta["sample_id"].values
    log.info(f"  Total samples: {len(meta)}")

    # ------------------------------------------------------------------
    # 2. Load global feature matrices
    # ------------------------------------------------------------------
    log.info("Loading global feature matrices...")
    ssgsea_df, tf_df, deconv_df = load_global_features(all_sample_ids)

    # ------------------------------------------------------------------
    # 3. Load LOCO NAS-VAE embeddings
    # ------------------------------------------------------------------
    log.info("Loading LOCO NAS-VAE embeddings...")
    loco_embeddings = {}
    for fold in NAS_DATASETS:
        emb_file = os.path.join(
            LOCO_EMB_DIR, f"nas_embeddings_loco_{fold}.csv"
        )
        if os.path.exists(emb_file):
            df = pd.read_csv(emb_file)
            df["sample_id"] = df["sample_id"].astype(str).str.strip()
            df = df.set_index("sample_id")
            loco_embeddings[fold] = df
            n_held = int(df["is_heldout"].sum()) if "is_heldout" in df.columns else "?"
            log.info(f"  {fold}: {df.shape[0]} samples, "
                     f"{n_held} held-out")
        else:
            log.warning(f"  Missing LOCO embeddings for {fold}: {emb_file}")

    # For fibrosis-only folds (GSE240729) that have no NAS labels and thus
    # no separate LOCO VAE model, we fall back to the full-model embeddings.
    # Load the full-model embeddings as fallback.
    full_emb_path = os.path.join(OUTDIR, "embeddings_all_samples.csv")
    full_emb_df = None
    if os.path.exists(full_emb_path):
        full_emb_df = pd.read_csv(full_emb_path)
        full_emb_df["sample_id"] = full_emb_df["sample_id"].astype(str).str.strip()
        full_emb_df = full_emb_df.set_index("sample_id")
        log.info(f"  Full-model embeddings fallback: {full_emb_df.shape}")

    # ------------------------------------------------------------------
    # 4. Define target configurations
    # ------------------------------------------------------------------
    # target_name -> (label_column, n_classes, cv_type)
    # cv_type: "nas" uses 5-fold NAS LOCO, "fib" uses 6-fold fibrosis LOCO
    target_cfg = {
        "nas_9class":  ("nas_score",  9, "nas"),
        "nas_4group":  ("nas_group4", 4, "nas"),
        "nas_ge5":     ("nas_ge5",    2, "nas"),
        "fib_stage":   ("fib_stage",  5, "fib"),
        "fib_ge3":     ("fib_ge3",    2, "fib"),
    }

    all_results = []

    # ------------------------------------------------------------------
    # 5. Run LOCO-CV for each target
    # ------------------------------------------------------------------
    for target_name, (label_col, n_classes, cv_type) in target_cfg.items():
        log.info(f"\n{'='*60}")
        log.info(f"Target: {target_name} ({label_col}, {n_classes} classes, "
                 f"{cv_type} LOCO)")
        log.info(f"{'='*60}")

        # Select fold assignments and datasets
        if cv_type == "nas":
            fold_col = "loco_fold_nas"
            fold_datasets = NAS_DATASETS
        else:
            fold_col = "loco_fold_fibrosis"
            fold_datasets = FIB_DATASETS

        # Filter to samples with valid labels and non-excluded folds
        eligible = meta[
            (meta[fold_col] != "excluded") & (meta[label_col] != MISSING)
        ].copy()
        log.info(f"  Eligible samples: {len(eligible)}")

        if len(eligible) < 50:
            log.warning(f"  Too few eligible samples ({len(eligible)}), skipping")
            continue

        # Check class distribution
        class_counts = eligible[label_col].value_counts().sort_index()
        log.info(f"  Class distribution: {class_counts.to_dict()}")

        for fold in fold_datasets:
            log.info(f"\n  --- Fold: {fold} ---")

            test_mask = eligible[fold_col] == fold
            train_mask = eligible[fold_col] != fold

            test_ids = eligible[test_mask]["sample_id"].values
            train_ids = eligible[train_mask]["sample_id"].values
            y_test = eligible[test_mask][label_col].values.astype(int)
            y_train = eligible[train_mask][label_col].values.astype(int)

            if len(test_ids) < 5:
                log.warning(f"    Skipping {fold}: only {len(test_ids)} test "
                            "samples")
                continue
            if len(train_ids) < 20:
                log.warning(f"    Skipping {fold}: only {len(train_ids)} train "
                            "samples")
                continue

            # Check that test fold has at least 2 classes for AUROC
            n_test_classes = len(np.unique(y_test))
            log.info(f"    Train: {len(train_ids)}, Test: {len(test_ids)} "
                     f"({n_test_classes} classes in test)")

            # Select embeddings: use LOCO-specific if available, else fallback
            if fold in loco_embeddings:
                emb_df = loco_embeddings[fold]
            elif full_emb_df is not None:
                log.info(f"    Using full-model embeddings fallback for {fold}")
                emb_df = full_emb_df
            else:
                log.warning(f"    No embeddings available for {fold}, skipping")
                continue

            # Assemble features
            try:
                X_train, X_test = assemble_fold_features(
                    fold_name=fold,
                    train_ids=train_ids,
                    test_ids=test_ids,
                    loco_emb_df=emb_df,
                    ssgsea_df=ssgsea_df,
                    tf_df=tf_df,
                    deconv_df=deconv_df,
                    meta_df=meta,
                )
            except Exception as e:
                log.error(f"    Feature assembly failed: {e}")
                continue

            # Handle edge case: only 1 class in training
            if len(np.unique(y_train)) < 2:
                log.warning(f"    Only 1 class in training set, skipping")
                continue

            # Train CORAL model
            try:
                result = train_coral_model(
                    X_train, y_train, X_test, y_test,
                    n_classes=n_classes,
                    target_name=target_name,
                    fold_name=fold,
                )
            except Exception as e:
                log.error(f"    Training failed: {e}")
                import traceback
                traceback.print_exc()
                continue

            log.info(f"    QWK={result['qwk']:.3f}  "
                     f"AUROC={result['auroc']:.3f}  "
                     f"MAE={result['mae']:.3f}")

            all_results.append({
                "target": target_name,
                "n_classes": n_classes,
                "cv_type": cv_type,
                "fold": fold,
                "n_train": len(train_ids),
                "n_test": len(test_ids),
                "qwk": result["qwk"],
                "auroc": result["auroc"],
                "mae": result["mae"],
            })

    # ------------------------------------------------------------------
    # 6. Save results
    # ------------------------------------------------------------------
    if not all_results:
        log.error("No results produced — check data availability")
        sys.exit(1)

    results_df = pd.DataFrame(all_results)
    results_df.to_csv(OUT_RESULTS, index=False)
    log.info(f"\nPer-fold results saved: {OUT_RESULTS}")
    log.info(f"  {len(results_df)} rows")

    # Summary: mean +/- std per target
    summary_rows = []
    for target_name in results_df["target"].unique():
        tdf = results_df[results_df["target"] == target_name]
        for metric in ["qwk", "auroc", "mae"]:
            vals = tdf[metric].dropna()
            if len(vals) > 0:
                summary_rows.append({
                    "target": target_name,
                    "n_classes": tdf["n_classes"].iloc[0],
                    "cv_type": tdf["cv_type"].iloc[0],
                    "metric": metric,
                    "mean": vals.mean(),
                    "std": vals.std(),
                    "min": vals.min(),
                    "max": vals.max(),
                    "n_folds": len(vals),
                })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT_SUMMARY, index=False)
    log.info(f"Summary saved: {OUT_SUMMARY}")

    # Print summary table
    print("\n" + "=" * 70)
    print("CORAL Multi-Scale Results Summary")
    print("=" * 70)
    for target_name in summary_df["target"].unique():
        tdf = summary_df[summary_df["target"] == target_name]
        print(f"\n  {target_name}:")
        for _, row in tdf.iterrows():
            print(f"    {row['metric']:>6s}: {row['mean']:.3f} +/- "
                  f"{row['std']:.3f}  [{row['min']:.3f} - {row['max']:.3f}]  "
                  f"(n={int(row['n_folds'])})")

    elapsed = time.time() - t0
    print(f"\nTotal time: {elapsed / 60:.1f} min")
    log.info(f"Done in {elapsed / 60:.1f} min")


if __name__ == "__main__":
    main()
