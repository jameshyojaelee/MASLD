#!/usr/bin/env python3
"""
92_multitask_semisupervised.py
Semi-supervised multi-task learning that leverages ALL resolution levels
simultaneously via a shared-backbone neural network with 7 task heads.

Architecture:
  Shared backbone:
    Linear(3000 -> 1024) -> BatchNorm -> ReLU -> Dropout(0.3)
    Linear(1024 ->  256) -> BatchNorm -> ReLU -> Dropout(0.3)
    Linear( 256 ->  128) -> BatchNorm -> ReLU -> Dropout(0.2)

  Task heads (each Linear(128 -> n_classes)):
    1. Disease binary      (128 -> 2) : 1,444 labels, weight 0.1
    2. NAS 3-class         (128 -> 3) :   660 labels, weight 0.3
    3. NAS 9-class ordinal (128 -> 9) :   660 labels, weight 0.5
    4. Steatosis ordinal   (128 -> 4) :    78 labels, weight 2.0
    5. Inflammation ordinal(128 -> 3) :    78 labels, weight 2.0
    6. Ballooning ordinal  (128 -> 3) :    78 labels, weight 2.0
    7. Fibrosis ordinal    (128 -> 5) : 1,128 labels, weight 0.3

Key features:
  - CORN ordinal encoding for tasks 3-7 (cumulative logit thresholds)
  - Label masking (loss only on samples with labels; -1 = missing)
  - Inverse-N task weighting
  - Hierarchical consistency regularisation (NAS 3-class vs NAS 9-class)

Training:
  - Adam lr=5e-4, batch=128, 300 epochs, early stopping patience 30
  - Cosine annealing LR scheduler
  - 5-fold NAS LOCO (primary) + 6-fold fibrosis LOCO (secondary)

Outputs (in results/staging_classifier/):
  - multitask_nas_results.csv      (per-fold, per-task metrics)
  - multitask_nas_summary.csv      (aggregated across folds)
  - multitask_training_curves.csv
  - multitask_model.pt             (best model weights from last fold)

Usage:  python 92_multitask_semisupervised.py
SLURM:  gpu, 1xL40S, 8 CPUs, 64G, 48h
Env:    micromamba activate rapids_singlecell
"""

import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
import h5py
from collections import defaultdict

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader

from sklearn.metrics import (
    roc_auc_score, accuracy_score, cohen_kappa_score
)

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)
torch.manual_seed(RANDOM_STATE)

# =============================================================================
# Paths
# =============================================================================
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
INT = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR = os.path.join(INT, "results")
OUTDIR = os.path.join(RDIR, "staging_classifier")
LOCO_GENE_DIR = os.path.join(OUTDIR, "loco_gene_lists")
os.makedirs(OUTDIR, exist_ok=True)

H5_PATH = os.path.join(OUTDIR, "prepared_data.h5")
META_PATH = os.path.join(OUTDIR, "modeling_metadata.csv")

# --- Device ---
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("=" * 70)
print("92: Semi-Supervised Multi-Task Learning (Upscaler)")
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

MISSING = -1

# LOCO fold definitions
NAS_DATASETS = ["GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066"]
FIB_DATASETS = ["GSE130970", "GSE135251", "GSE162694", "GSE174478",
                "GSE193066", "GSE240729"]

# Hierarchical consistency: NAS 3-class ↔ NAS 9-class mapping
# nas_3class: 0=Low (NAS 0-2), 1=Mid (NAS 3-4), 2=High (NAS 5-8)
NAS9_TO_NAS3 = {0: 0, 1: 0, 2: 0, 3: 1, 4: 1, 5: 2, 6: 2, 7: 2, 8: 2}

# =============================================================================
# Load data
# =============================================================================
print("Loading data...")

if not os.path.exists(H5_PATH):
    sys.exit(f"ERROR: {H5_PATH} not found. Run Script 62 first.")
if not os.path.exists(META_PATH):
    sys.exit(f"ERROR: {META_PATH} not found. Run Script 89 first.")

with h5py.File(H5_PATH, "r") as h5:
    X_all = h5["rank_expression"][:].T.astype(np.float32)  # (samples, genes)
    gene_names = np.array([g.decode() if isinstance(g, bytes) else g
                           for g in h5["gene_names"][:]])
    sample_ids = np.array([s.decode() if isinstance(s, bytes) else s
                           for s in h5["sample_ids"][:]])

N_SAMPLES, N_GENES = X_all.shape
print(f"  Expression: {N_SAMPLES} samples x {N_GENES} genes")

# Build gene name -> index mapping for fold-specific subsetting
_gene_name_to_idx_92 = {g: i for i, g in enumerate(gene_names)}

# Load global feature candidates as fallback
_global_fc_path_92 = os.path.join(OUTDIR, "feature_candidates_3000.csv")
if os.path.exists(_global_fc_path_92):
    _global_fc_92 = pd.read_csv(_global_fc_path_92)
    _feature_genes_global_92 = _global_fc_92["gene"].tolist()
else:
    _feature_genes_global_92 = list(gene_names)


def load_fold_gene_indices_92(fold_id):
    """Return column indices into X_all for the per-fold gene list.

    Falls back to the global feature_candidates_3000.csv if the fold-specific
    file does not exist.
    """
    fold_gene_file = os.path.join(
        LOCO_GENE_DIR, f"fold_{fold_id}_feature_candidates_3000.csv")
    if os.path.exists(fold_gene_file):
        fold_genes = pd.read_csv(fold_gene_file)["gene"].tolist()
        print(f"    Using per-fold gene list for fold: {fold_id}")
    else:
        fold_genes = _feature_genes_global_92
        print(f"    Fold gene list not found for {fold_id} -- using global fallback")
    return np.array([_gene_name_to_idx_92[g] for g in fold_genes
                     if g in _gene_name_to_idx_92], dtype=np.int64)

# Load updated metadata (with steatosis, inflammation, ballooning, nas_3class)
meta_df = pd.read_csv(META_PATH)
# Ensure sample_id is string in both
meta_df["sample_id"] = meta_df["sample_id"].astype(str).str.strip()
sample_ids_str = [str(s).strip() for s in sample_ids]
# Align to H5 sample order
meta_df = meta_df.set_index("sample_id")
matched_ids = [s for s in sample_ids_str if s in meta_df.index]
print(f"  Matched {len(matched_ids)} of {len(sample_ids_str)} samples to metadata")
if len(matched_ids) < len(sample_ids_str):
    # Try without leading/trailing whitespace or version differences
    meta_set = set(meta_df.index)
    h5_set = set(sample_ids_str)
    print(f"  Meta samples: {len(meta_set)}, H5 samples: {len(h5_set)}")
    print(f"  Meta sample example: '{list(meta_set)[:3]}'")
    print(f"  H5 sample example: '{sample_ids_str[:3]}'")
meta_df = meta_df.loc[matched_ids].reset_index()
assert len(meta_df) > 0, f"No metadata matches! Check sample ID format."
if len(meta_df) < N_SAMPLES:
    print(f"  WARNING: Only {len(meta_df)}/{N_SAMPLES} samples matched metadata")
    # Subset expression to matched samples
    match_idx = [i for i, s in enumerate(sample_ids_str) if s in set(matched_ids)]
    X_all = X_all[match_idx]
    sample_ids = np.array(matched_ids)
    N_SAMPLES = len(meta_df)
    print(f"  Subsetted expression to {N_SAMPLES} samples")

# =============================================================================
# Build label arrays  (each shape (N_SAMPLES,), dtype int64, -1 = missing)
# =============================================================================
print("\nBuilding label arrays...")


def safe_int_label(series, valid_range=None):
    """Convert to int64 array; anything outside valid_range -> MISSING."""
    arr = series.fillna(MISSING).values.astype(np.int64)
    if valid_range is not None:
        lo, hi = valid_range
        arr[(arr < lo) | (arr > hi)] = MISSING
    return arr


labels = {}

# 1. Disease binary (0=Control, 1=Disease) -- available for all 1,444 samples
labels["disease"] = safe_int_label(meta_df["is_disease"], (0, 1))

# 2. NAS 3-class (0=Low, 1=Mid, 2=High)
labels["nas_3class"] = safe_int_label(meta_df["nas_3class"], (0, 2))

# 3. NAS 9-class ordinal (0-8); cap raw NAS at 8
nas_raw = meta_df["nas_score"].fillna(-1).values.astype(np.int64)
nas_raw[nas_raw > 8] = 8
nas_raw[nas_raw < 0] = MISSING
labels["nas_9class"] = nas_raw

# 4. Steatosis ordinal (0-3)
labels["steatosis"] = safe_int_label(meta_df["steatosis_grade"], (0, 3))

# 5. Inflammation ordinal (0-2)
labels["inflammation"] = safe_int_label(
    meta_df["lobular_inflammation_grade"], (0, 2))

# 6. Ballooning ordinal (0-2)
labels["ballooning"] = safe_int_label(meta_df["ballooning_grade"], (0, 2))

# 7. Fibrosis ordinal (0-4)
labels["fibrosis"] = safe_int_label(meta_df["fib_stage"], (0, 4))

print("Label availability:")
for t, arr in labels.items():
    n_valid = int(np.sum(arr >= 0))
    n_classes = TASK_CFG[t][0]
    dist = dict(zip(*np.unique(arr[arr >= 0], return_counts=True))) if n_valid > 0 else {}
    print(f"  {t:20s}: {n_valid:5d} labels, {n_classes} classes, dist={dist}")

datasets = meta_df["dataset"].values


# =============================================================================
# CORN ordinal encoding utilities
# =============================================================================
class CORNLayer(nn.Module):
    """Conditional Ordinal Regression Network (CORN) head.

    For K ordinal classes, predicts K-1 cumulative logits.
    P(Y > k | Y >= k) is modelled independently for each threshold k.
    Class probabilities are derived via the chain rule:
        P(Y = 0) = 1 - sigma(z_0)
        P(Y = k) = prod_{j<k} sigma(z_j) * (1 - sigma(z_k))  for 0 < k < K-1
        P(Y = K-1) = prod_{j<K-1} sigma(z_j)
    """

    def __init__(self, in_features: int, n_classes: int):
        super().__init__()
        self.n_classes = n_classes
        # One logit per threshold (K-1 thresholds for K classes)
        self.fc = nn.Linear(in_features, n_classes - 1)

    def forward(self, x):
        """Return raw logits of shape (batch, n_classes - 1)."""
        return self.fc(x)


def corn_logits_to_probs(logits: torch.Tensor, n_classes: int) -> torch.Tensor:
    """Convert CORN logits (batch, K-1) -> class probabilities (batch, K).

    Uses the conditional chain-rule formulation:
        sigma_k = sigmoid(logit_k) = P(Y > k | Y >= k)
        P(Y = 0) = 1 - sigma_0
        P(Y = k) = (prod_{j<k} sigma_j) * (1 - sigma_k)
        P(Y = K-1) = prod_{j<K-1} sigma_j
    """
    sigmas = torch.sigmoid(logits)  # (batch, K-1)
    batch = sigmas.shape[0]
    probs = torch.zeros(batch, n_classes, device=logits.device)

    # P(Y=0) = 1 - sigma_0
    probs[:, 0] = 1.0 - sigmas[:, 0]

    # P(Y=k) for 1 <= k <= K-2
    cum_prod = sigmas[:, 0].clone()
    for k in range(1, n_classes - 1):
        probs[:, k] = cum_prod * (1.0 - sigmas[:, k])
        cum_prod = cum_prod * sigmas[:, k]

    # P(Y=K-1) = prod of all sigmas
    probs[:, n_classes - 1] = cum_prod

    # Numerical safety: clamp and re-normalise
    probs = probs.clamp(min=1e-7)
    probs = probs / probs.sum(dim=1, keepdim=True)
    return probs


def corn_loss(logits: torch.Tensor, targets: torch.Tensor,
              n_classes: int) -> torch.Tensor:
    """Compute CORN loss (masked, ignoring targets == -1).

    For each threshold k, the loss is binary cross-entropy for the event
    "Y > k" conditioned on samples where Y >= k.
    """
    mask = targets >= 0
    if mask.sum() == 0:
        return torch.tensor(0.0, device=logits.device, requires_grad=True)

    logits_v = logits[mask]   # (n_valid, K-1)
    targets_v = targets[mask]  # (n_valid,)
    loss = torch.tensor(0.0, device=logits.device)
    n_thresholds = 0

    for k in range(n_classes - 1):
        # Samples with Y >= k
        cond_mask = targets_v >= k
        if cond_mask.sum() < 2:
            continue
        # Binary target: 1 if Y > k, else 0
        binary_target = (targets_v[cond_mask] > k).float()
        bce = nn.functional.binary_cross_entropy_with_logits(
            logits_v[cond_mask, k], binary_target, reduction="mean"
        )
        loss = loss + bce
        n_thresholds += 1

    if n_thresholds > 0:
        loss = loss / n_thresholds
    return loss


# =============================================================================
# Masked cross-entropy for non-ordinal (disease) task
# =============================================================================
def masked_cross_entropy(logits, targets, weight=None):
    """Standard cross-entropy, ignoring targets == -1."""
    mask = targets >= 0
    if mask.sum() == 0:
        return torch.tensor(0.0, device=logits.device, requires_grad=True)
    return nn.functional.cross_entropy(
        logits[mask], targets[mask], weight=weight)


# =============================================================================
# Inverse-frequency class weights
# =============================================================================
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
# QWK helper
# =============================================================================
def quadratic_weighted_kappa(y_true, y_pred):
    """Manual QWK implementation."""
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


# =============================================================================
# Dataset
# =============================================================================
class SemiSupervisedDataset(Dataset):
    """Returns expression features + per-task label dict."""

    def __init__(self, X, label_dict):
        self.X = torch.tensor(X, dtype=torch.float32)
        self.labels = {
            k: torch.tensor(v, dtype=torch.long)
            for k, v in label_dict.items()
        }

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        return self.X[idx], {k: v[idx] for k, v in self.labels.items()}


# =============================================================================
# Multi-Task Model
# =============================================================================
class MultiTaskSemiSupervised(nn.Module):
    """Shared backbone with BatchNorm + 7 task heads.

    Non-ordinal heads output class logits; ordinal heads use CORN.
    """

    def __init__(self, n_features: int):
        super().__init__()

        self.backbone = nn.Sequential(
            nn.Linear(n_features, 1024),
            nn.BatchNorm1d(1024),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(1024, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(256, 128),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.Dropout(0.2),
        )

        # Non-ordinal head
        self.head_disease = nn.Linear(128, 2)

        # Ordinal (CORN) heads
        self.head_nas_3class = CORNLayer(128, 3)
        self.head_nas_9class = CORNLayer(128, 9)
        self.head_steatosis = CORNLayer(128, 4)
        self.head_inflammation = CORNLayer(128, 3)
        self.head_ballooning = CORNLayer(128, 3)
        self.head_fibrosis = CORNLayer(128, 5)

    def forward(self, x):
        z = self.backbone(x)
        return {
            "disease": self.head_disease(z),
            "nas_3class": self.head_nas_3class(z),
            "nas_9class": self.head_nas_9class(z),
            "steatosis": self.head_steatosis(z),
            "inflammation": self.head_inflammation(z),
            "ballooning": self.head_ballooning(z),
            "fibrosis": self.head_fibrosis(z),
        }


# =============================================================================
# Hierarchical consistency regularisation
# =============================================================================
def hierarchical_consistency_loss(logits_3class, logits_9class):
    """Penalise when NAS 3-class and 9-class predictions disagree.

    Marginalise NAS 9-class probabilities down to 3-class and compute
    KL divergence against the NAS 3-class predicted distribution.
    """
    # Convert CORN logits to probabilities
    probs_3 = corn_logits_to_probs(logits_3class, 3)  # (batch, 3)
    probs_9 = corn_logits_to_probs(logits_9class, 9)  # (batch, 9)

    # Marginalise 9-class -> 3-class
    # 0-2 -> class 0, 3-4 -> class 1, 5-8 -> class 2
    probs_9_marg = torch.zeros_like(probs_3)
    probs_9_marg[:, 0] = probs_9[:, 0:3].sum(dim=1)
    probs_9_marg[:, 1] = probs_9[:, 3:5].sum(dim=1)
    probs_9_marg[:, 2] = probs_9[:, 5:9].sum(dim=1)

    # Clamp for numerical stability
    probs_3 = probs_3.clamp(min=1e-7)
    probs_9_marg = probs_9_marg.clamp(min=1e-7)

    # Symmetric KL: 0.5 * (KL(p3 || p9m) + KL(p9m || p3))
    kl_fwd = (probs_3 * (probs_3.log() - probs_9_marg.log())).sum(dim=1)
    kl_rev = (probs_9_marg * (probs_9_marg.log() - probs_3.log())).sum(dim=1)
    return 0.5 * (kl_fwd + kl_rev).mean()


# =============================================================================
# Training
# =============================================================================
def train_one_fold(X_train, lab_train, X_val, lab_val,
                   n_epochs=300, lr=5e-4, batch_size=128, patience=30,
                   consistency_weight=0.1, fold_name="fold"):
    """Train multi-task model for one LOCO fold; return predictions + curves."""

    train_ds = SemiSupervisedDataset(X_train, lab_train)
    val_ds = SemiSupervisedDataset(X_val, lab_val)

    train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True,
                              drop_last=False, num_workers=0, pin_memory=True)
    val_loader = DataLoader(val_ds, batch_size=batch_size, shuffle=False,
                            num_workers=0, pin_memory=True)

    # Class weights for the non-ordinal disease head
    cw_disease = compute_class_weights(lab_train["disease"], 2, device)

    model = MultiTaskSemiSupervised(X_train.shape[1]).to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=n_epochs, eta_min=1e-6)

    best_val_loss = float("inf")
    best_state = None
    patience_counter = 0
    history = []

    for epoch in range(n_epochs):
        # ---- Train ----
        model.train()
        train_losses = defaultdict(float)
        n_batches = 0

        for X_batch, y_batch in train_loader:
            X_batch = X_batch.to(device, non_blocking=True)
            y_batch = {k: v.to(device, non_blocking=True)
                       for k, v in y_batch.items()}

            logits = model(X_batch)
            loss = torch.tensor(0.0, device=device)

            # Disease: standard CE
            task_loss = masked_cross_entropy(
                logits["disease"], y_batch["disease"], weight=cw_disease)
            loss = loss + TASK_CFG["disease"][1] * task_loss
            train_losses["disease"] += task_loss.item()

            # Ordinal tasks: CORN loss
            for task in ["nas_3class", "nas_9class", "steatosis",
                         "inflammation", "ballooning", "fibrosis"]:
                n_cls, tw, _ = TASK_CFG[task]
                tl = corn_loss(logits[task], y_batch[task], n_cls)
                loss = loss + tw * tl
                train_losses[task] += tl.item()

            # Hierarchical consistency (on ALL samples, not just labelled)
            hc_loss = hierarchical_consistency_loss(
                logits["nas_3class"], logits["nas_9class"])
            loss = loss + consistency_weight * hc_loss
            train_losses["consistency"] += hc_loss.item()

            optimizer.zero_grad()
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()
            n_batches += 1

        scheduler.step()

        # ---- Validate ----
        model.eval()
        val_losses = defaultdict(float)
        n_val_batches = 0
        with torch.no_grad():
            for X_batch, y_batch in val_loader:
                X_batch = X_batch.to(device, non_blocking=True)
                y_batch = {k: v.to(device, non_blocking=True)
                           for k, v in y_batch.items()}
                logits = model(X_batch)

                # Disease
                vl = masked_cross_entropy(
                    logits["disease"], y_batch["disease"], weight=cw_disease)
                val_losses["disease"] += vl.item()

                # Ordinal tasks
                for task in ["nas_3class", "nas_9class", "steatosis",
                             "inflammation", "ballooning", "fibrosis"]:
                    n_cls = TASK_CFG[task][0]
                    vl = corn_loss(logits[task], y_batch[task], n_cls)
                    val_losses[task] += vl.item()

                # Consistency
                hc = hierarchical_consistency_loss(
                    logits["nas_3class"], logits["nas_9class"])
                val_losses["consistency"] += hc.item()
                n_val_batches += 1

        # Weighted total validation loss
        total_val = sum(
            TASK_CFG[t][1] * val_losses[t] / max(n_val_batches, 1)
            for t in TASK_CFG
        ) + consistency_weight * val_losses["consistency"] / max(n_val_batches, 1)

        row = {
            "fold": fold_name,
            "epoch": epoch + 1,
            "lr": scheduler.get_last_lr()[0],
            "total_val_loss": total_val,
        }
        for t in list(TASK_CFG.keys()) + ["consistency"]:
            row[f"train_{t}"] = train_losses[t] / max(n_batches, 1)
            row[f"val_{t}"] = val_losses[t] / max(n_val_batches, 1)
        history.append(row)

        # Early stopping
        if total_val < best_val_loss:
            best_val_loss = total_val
            best_state = {k: v.cpu().clone()
                          for k, v in model.state_dict().items()}
            patience_counter = 0
        else:
            patience_counter += 1

        if patience_counter >= patience:
            print(f"  [{fold_name}] Early stopping at epoch {epoch + 1}")
            break

        if (epoch + 1) % 25 == 0 or epoch == 0:
            msg_parts = [f"e={epoch+1:3d}"]
            for t in TASK_CFG:
                msg_parts.append(
                    f"{t[:6]}={train_losses[t]/max(n_batches,1):.3f}")
            msg_parts.append(f"val={total_val:.4f}")
            print(f"  [{fold_name}] {' | '.join(msg_parts)}")

    # Load best model
    if best_state is not None:
        model.load_state_dict(best_state)
    model.to(device)

    # ---- Predict on validation set ----
    model.eval()
    all_logits = {t: [] for t in TASK_CFG}
    with torch.no_grad():
        for X_batch, _ in val_loader:
            X_batch = X_batch.to(device, non_blocking=True)
            out = model(X_batch)
            for t in TASK_CFG:
                all_logits[t].append(out[t].cpu())

    predictions = {}
    for t in TASK_CFG:
        cat = torch.cat(all_logits[t], dim=0)
        n_cls, _, is_ord = TASK_CFG[t]
        if is_ord:
            probs = corn_logits_to_probs(cat, n_cls).numpy()
        else:
            probs = torch.softmax(cat, dim=1).numpy()
        predictions[t] = probs

    return predictions, pd.DataFrame(history), model


# =============================================================================
# Evaluation helpers
# =============================================================================
def eval_ordinal(y_true, probs, task_name, n_classes):
    """Evaluate an ordinal task: QWK, MAE, accuracy."""
    valid = y_true >= 0
    n_valid = int(valid.sum())
    if n_valid < 5:
        return {"task": task_name, "n_valid": n_valid,
                "qwk": np.nan, "mae": np.nan, "acc": np.nan,
                "adj_acc": np.nan}

    y_v = y_true[valid].astype(int)
    p_v = probs[valid]
    y_pred = p_v.argmax(axis=1)

    qwk = quadratic_weighted_kappa(y_v, y_pred)
    mae = float(np.mean(np.abs(y_v - y_pred)))
    acc = float(accuracy_score(y_v, y_pred))
    adj_acc = float(np.mean(np.abs(y_v - y_pred) <= 1))

    return {"task": task_name, "n_valid": n_valid,
            "qwk": qwk, "mae": mae, "acc": acc, "adj_acc": adj_acc}


def eval_binary_auroc(y_true, probs, task_name, threshold_col, threshold_val):
    """Binarise ordinal labels at threshold and compute AUROC."""
    valid = y_true >= 0
    n_valid = int(valid.sum())
    if n_valid < 10:
        return {"task": task_name, "metric": f"{task_name}_auroc",
                "n_valid": n_valid, "value": np.nan}
    y_v = y_true[valid]
    p_v = probs[valid]
    # P(Y >= threshold) = sum of probs for classes >= threshold
    prob_ge = p_v[:, threshold_val:].sum(axis=1)
    binary = (y_v >= threshold_val).astype(int)
    if len(np.unique(binary)) < 2:
        return {"task": task_name, "metric": f"{task_name}_auroc",
                "n_valid": n_valid, "value": np.nan}
    auroc = roc_auc_score(binary, prob_ge)
    return {"task": task_name, "metric": f"{task_name}_auroc",
            "n_valid": n_valid, "value": auroc}


def eval_fold(predictions, lab_val, fold_name):
    """Compute all evaluation metrics for one LOCO fold."""
    rows = []

    # --- NAS 9-class ---
    r = eval_ordinal(lab_val["nas_9class"], predictions["nas_9class"],
                     "nas_9class", 9)
    r["fold"] = fold_name
    rows.append(r)

    # NAS 4-group (0-2=Low, 3-4=Mid, 5-6=High, 7-8=VHigh)
    valid_9 = lab_val["nas_9class"] >= 0
    if valid_9.sum() >= 5:
        y9 = lab_val["nas_9class"][valid_9].astype(int)
        p9 = predictions["nas_9class"][valid_9]
        # Marginalise to 4 groups
        y4 = np.zeros_like(y9)
        y4[y9 <= 2] = 0
        y4[(y9 >= 3) & (y9 <= 4)] = 1
        y4[(y9 >= 5) & (y9 <= 6)] = 2
        y4[y9 >= 7] = 3

        p4 = np.zeros((len(p9), 4))
        p4[:, 0] = p9[:, 0:3].sum(axis=1)
        p4[:, 1] = p9[:, 3:5].sum(axis=1)
        p4[:, 2] = p9[:, 5:7].sum(axis=1)
        p4[:, 3] = p9[:, 7:9].sum(axis=1)
        yp4 = p4.argmax(axis=1)
        qwk4 = quadratic_weighted_kappa(y4, yp4)
        rows.append({"task": "nas_4group", "fold": fold_name,
                      "n_valid": int(valid_9.sum()), "qwk": qwk4,
                      "mae": np.nan, "acc": float(accuracy_score(y4, yp4)),
                      "adj_acc": np.nan})

    # --- NAS 3-class ---
    r = eval_ordinal(lab_val["nas_3class"], predictions["nas_3class"],
                     "nas_3class", 3)
    r["fold"] = fold_name
    rows.append(r)

    # --- NAS >= 5 AUROC (from 9-class) ---
    r_auroc = eval_binary_auroc(
        lab_val["nas_9class"], predictions["nas_9class"],
        "nas_ge5", "nas_9class", 5)
    r_auroc["fold"] = fold_name
    rows.append({"task": "nas_ge5_auroc", "fold": fold_name,
                  "n_valid": r_auroc["n_valid"],
                  "qwk": np.nan, "mae": np.nan,
                  "acc": np.nan, "adj_acc": np.nan,
                  "auroc": r_auroc["value"]})

    # --- Steatosis (only GSE130970 samples) ---
    r = eval_ordinal(lab_val["steatosis"], predictions["steatosis"],
                     "steatosis", 4)
    r["fold"] = fold_name
    rows.append(r)

    # --- Inflammation ---
    r = eval_ordinal(lab_val["inflammation"], predictions["inflammation"],
                     "inflammation", 3)
    r["fold"] = fold_name
    rows.append(r)

    # --- Ballooning ---
    r = eval_ordinal(lab_val["ballooning"], predictions["ballooning"],
                     "ballooning", 3)
    r["fold"] = fold_name
    rows.append(r)

    # --- Fibrosis ---
    r = eval_ordinal(lab_val["fibrosis"], predictions["fibrosis"],
                     "fibrosis", 5)
    r["fold"] = fold_name
    rows.append(r)

    # --- Disease ---
    valid_d = lab_val["disease"] >= 0
    if valid_d.sum() >= 5:
        yd = lab_val["disease"][valid_d].astype(int)
        pd_ = predictions["disease"][valid_d]
        ypd = pd_.argmax(axis=1)
        acc_d = float(accuracy_score(yd, ypd))
        try:
            auroc_d = roc_auc_score(yd, pd_[:, 1]) if len(np.unique(yd)) > 1 else np.nan
        except ValueError:
            auroc_d = np.nan
        rows.append({"task": "disease", "fold": fold_name,
                      "n_valid": int(valid_d.sum()),
                      "qwk": np.nan, "mae": np.nan,
                      "acc": acc_d, "adj_acc": np.nan,
                      "auroc": auroc_d})

    return rows


# =============================================================================
# Main LOCO evaluation
# =============================================================================
def run_loco(fold_type, fold_datasets, all_X, all_labels, all_datasets):
    """Run LOCO cross-validation.

    Args:
        fold_type: "NAS" or "FIB"
        fold_datasets: list of dataset names forming folds
        all_X: (N, G) expression array
        all_labels: dict of label arrays
        all_datasets: array of dataset names per sample

    Returns:
        results_rows: list of metric dicts
        curve_df: training curves across all folds
        last_model: model from the final fold (for saving)
    """
    print(f"\n{'='*60}")
    print(f"LOCO Evaluation: {fold_type} ({len(fold_datasets)} folds)")
    print(f"{'='*60}")

    all_results = []
    all_curves = []
    last_model = None

    for fold_idx, held_out in enumerate(fold_datasets):
        print(f"\n--- Fold {fold_idx+1}/{len(fold_datasets)}: "
              f"held-out = {held_out} ---")

        # Load per-fold gene list and subset features (removes leakage)
        fold_gene_idx = load_fold_gene_indices_92(held_out)
        X_fold = all_X[:, fold_gene_idx]
        n_features_fold = X_fold.shape[1]
        print(f"  Features for this fold: {n_features_fold}")

        # Split: train = all samples NOT in held-out dataset
        # (use ALL samples for training, not just those in fold_datasets,
        #  because semi-supervised leverages unlabelled samples too)
        val_mask = all_datasets == held_out
        train_mask = ~val_mask

        X_train = X_fold[train_mask]
        X_val = X_fold[val_mask]

        lab_train = {t: arr[train_mask] for t, arr in all_labels.items()}
        lab_val = {t: arr[val_mask] for t, arr in all_labels.items()}

        n_train = int(train_mask.sum())
        n_val = int(val_mask.sum())
        print(f"  Train: {n_train} samples, Val: {n_val} samples")

        # Print per-task label counts in train/val
        for t in TASK_CFG:
            nt = int(np.sum(lab_train[t] >= 0))
            nv = int(np.sum(lab_val[t] >= 0))
            print(f"    {t:20s}: train={nt:5d}, val={nv:4d}")

        fold_name = f"{fold_type}_{held_out}"

        predictions, curves, model = train_one_fold(
            X_train, lab_train, X_val, lab_val, fold_name=fold_name)
        last_model = model

        fold_results = eval_fold(predictions, lab_val, fold_name)
        all_results.extend(fold_results)
        all_curves.append(curves)

        # Print key metrics for this fold
        for r in fold_results:
            parts = [f"{r['task']:20s}"]
            if not np.isnan(r.get("qwk", np.nan)):
                parts.append(f"QWK={r['qwk']:.3f}")
            if not np.isnan(r.get("mae", np.nan)):
                parts.append(f"MAE={r['mae']:.3f}")
            if not np.isnan(r.get("acc", np.nan)):
                parts.append(f"Acc={r['acc']:.3f}")
            if not np.isnan(r.get("auroc", np.nan)):
                parts.append(f"AUROC={r['auroc']:.3f}")
            print(f"    {'  '.join(parts)}")

    curve_df = pd.concat(all_curves, ignore_index=True)
    return all_results, curve_df, last_model


# =============================================================================
# Main
# =============================================================================
def main():
    t_start = time.time()

    # ---- NAS LOCO (primary) ----
    nas_results, nas_curves, nas_model = run_loco(
        "NAS", NAS_DATASETS, X_all, labels, datasets)

    # ---- Fibrosis LOCO (secondary) ----
    fib_results, fib_curves, fib_model = run_loco(
        "FIB", FIB_DATASETS, X_all, labels, datasets)

    # ---- Combine and save results ----
    all_results = nas_results + fib_results
    results_df = pd.DataFrame(all_results)
    results_path = os.path.join(OUTDIR, "multitask_nas_results.csv")
    results_df.to_csv(results_path, index=False)
    print(f"\nPer-fold results saved: {results_path}")

    # ---- Aggregate summary ----
    print(f"\n{'='*70}")
    print("AGGREGATED RESULTS")
    print(f"{'='*70}")

    summary_rows = []
    for fold_type in ["NAS", "FIB"]:
        prefix = f"{fold_type}_"
        sub = results_df[results_df["fold"].str.startswith(prefix)]
        for task in sub["task"].unique():
            tsub = sub[sub["task"] == task]
            row = {"fold_type": fold_type, "task": task,
                   "n_folds": len(tsub)}
            for metric in ["qwk", "mae", "acc", "adj_acc", "auroc"]:
                if metric in tsub.columns:
                    vals = tsub[metric].dropna()
                    if len(vals) > 0:
                        row[f"{metric}_mean"] = vals.mean()
                        row[f"{metric}_std"] = vals.std()
                    else:
                        row[f"{metric}_mean"] = np.nan
                        row[f"{metric}_std"] = np.nan
            summary_rows.append(row)

            # Print
            parts = [f"[{fold_type}] {task:20s}"]
            for m in ["qwk", "mae", "acc", "auroc"]:
                key = f"{m}_mean"
                if key in row and not np.isnan(row.get(key, np.nan)):
                    parts.append(
                        f"{m}={row[key]:.3f}+/-{row.get(f'{m}_std',0):.3f}")
            print(f"  {'  '.join(parts)}")

    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(OUTDIR, "multitask_nas_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"\nSummary saved: {summary_path}")

    # ---- Training curves ----
    curves_df = pd.concat([nas_curves, fib_curves], ignore_index=True)
    curves_path = os.path.join(OUTDIR, "multitask_training_curves.csv")
    curves_df.to_csv(curves_path, index=False)
    print(f"Training curves saved: {curves_path}")

    # ---- Save best model (from last NAS fold) ----
    model_path = os.path.join(OUTDIR, "multitask_model.pt")
    if nas_model is not None:
        torch.save(nas_model.state_dict(), model_path)
        print(f"Model weights saved: {model_path}")

    elapsed = time.time() - t_start
    print(f"\nTotal time: {elapsed/60:.1f} minutes")
    print(f"\n=== 92_multitask_semisupervised.py completed ===")


if __name__ == "__main__":
    main()
