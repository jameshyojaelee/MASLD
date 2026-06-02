#!/usr/bin/env python3
"""
72_multimodal_hooks.py
Multi-modal architecture with fusion demo.

Defines modality-specific encoders and a fusion model. Demonstrates
bulk RNA + atlas feature fusion vs bulk-only for fibrosis prediction.

Architecture:
  BulkRNAEncoder:      Pre-trained VAE encoder (frozen) -> 64-dim
  AtlasFeatureEncoder:  Linear(N->64)->ReLU->Dropout->Linear(64->32) -> 32-dim
  FusionModel:          Concat(64+32=96)->Linear(96->64)->ReLU->Dropout->Linear(64->n_classes)

  Placeholder encoders: ScRNAEncoder, ATACEncoder, ProteomicsEncoder (for future use)

Demo experiment:
  Load pre-trained VAE encoder, freeze weights
  Load atlas features from HDF5
  Train fusion: bulk RNA embedding + atlas features -> fibrosis prediction
  Compare: bulk-only vs fusion via LOCO validation

Input:
  - results/staging_classifier/prepared_data.h5 (expression + atlas features)
  - results/staging_classifier/vae_model_full.pt (trained VAE encoder)

Output: results/staging_classifier/
  - multimodal_architecture.py (importable module with all encoder classes)
  - fusion_demo_results.csv (bulk-only vs fusion AUROC per fold)
  - modality_contribution.csv (ablation: which modality contributes more)

SLURM: gpu partition, 1xL40S, 4 CPUs, 32GB RAM, 48h
Env:    micromamba activate rapids_singlecell
"""

import os
import sys
import time
import warnings
import copy
import numpy as np
import pandas as pd
import h5py
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import roc_auc_score, balanced_accuracy_score, cohen_kappa_score
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
INT = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
OUTDIR = os.path.join(INT, "results/staging_classifier")
os.makedirs(OUTDIR, exist_ok=True)

H5_PATH = os.path.join(OUTDIR, "prepared_data.h5")
VAE_PATH = os.path.join(OUTDIR, "vae_model_full.pt")
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

print("=" * 60)
print("72: Multi-Modal Fusion Architecture + Demo")
print("=" * 60)
print(f"Device: {DEVICE}")
print()


# ===========================================================================
# ARCHITECTURE DEFINITIONS
# ===========================================================================
# These are also saved as a standalone importable module (multimodal_architecture.py)

class MultiModalEncoder(nn.Module):
    """Base class for modality-specific encoders."""

    def __init__(self, output_dim: int):
        super().__init__()
        self.output_dim = output_dim

    def encode(self, x):
        raise NotImplementedError

    def forward(self, x):
        return self.encode(x)


class BulkRNAEncoder(MultiModalEncoder):
    """
    Uses pre-trained VAE encoder to map 3000 genes -> 64-dim embedding.
    VAE encoder weights can be frozen or fine-tuned.
    """

    def __init__(self, n_genes=3000, latent_dim=64, dropout=0.3, freeze=True):
        super().__init__(output_dim=latent_dim)
        self.encoder = nn.Sequential(
            nn.Linear(n_genes, 2048),
            nn.BatchNorm1d(2048),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(2048, 512),
            nn.BatchNorm1d(512),
            nn.ReLU(),
            nn.Dropout(dropout),
        )
        self.fc_mu = nn.Linear(512, latent_dim)
        self.freeze = freeze

    def load_from_vae(self, vae_state_dict):
        """Load encoder weights from a trained ConditionalVAE state dict."""
        # Map VAE keys to our keys
        mapping = {}
        for key in vae_state_dict:
            if key.startswith("encoder."):
                mapping[key] = vae_state_dict[key]
            elif key == "fc_mu.weight":
                mapping[key] = vae_state_dict[key]
            elif key == "fc_mu.bias":
                mapping[key] = vae_state_dict[key]

        self.load_state_dict(mapping, strict=False)

        if self.freeze:
            for param in self.parameters():
                param.requires_grad = False

    def encode(self, x):
        h = self.encoder(x)
        return self.fc_mu(h)


class AtlasFeatureEncoder(MultiModalEncoder):
    """
    Maps N atlas numeric columns -> 32-dim embedding.
    Linear(N->64)->ReLU->Dropout->Linear(64->32)
    """

    def __init__(self, n_features, output_dim=32, dropout=0.2):
        super().__init__(output_dim=output_dim)
        self.net = nn.Sequential(
            nn.Linear(n_features, 64),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(64, output_dim),
        )

    def encode(self, x):
        return self.net(x)


class ScRNAEncoder(MultiModalEncoder):
    """Placeholder: aggregated cell-type proportion/expression embeddings."""

    def __init__(self, n_celltypes=10, output_dim=32):
        super().__init__(output_dim=output_dim)
        self.net = nn.Sequential(
            nn.Linear(n_celltypes, 64),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(64, output_dim),
        )

    def encode(self, x):
        return self.net(x)


class ATACEncoder(MultiModalEncoder):
    """Placeholder: chromatin accessibility peak features."""

    def __init__(self, n_peaks=500, output_dim=32):
        super().__init__(output_dim=output_dim)
        self.net = nn.Sequential(
            nn.Linear(n_peaks, 256),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(256, output_dim),
        )

    def encode(self, x):
        return self.net(x)


class ProteomicsEncoder(MultiModalEncoder):
    """Placeholder: protein abundance features."""

    def __init__(self, n_proteins=200, output_dim=32):
        super().__init__(output_dim=output_dim)
        self.net = nn.Sequential(
            nn.Linear(n_proteins, 128),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(128, output_dim),
        )

    def encode(self, x):
        return self.net(x)


class FusionModel(nn.Module):
    """
    Concatenates modality embeddings -> shared MLP -> prediction heads.

    Takes a list of (encoder, input_tensor) pairs.
    The encoders produce embeddings that are concatenated and fed through
    a shared classification MLP.
    """

    def __init__(self, encoder_dims, n_classes, hidden_dim=64, dropout=0.3):
        """
        Args:
            encoder_dims: list of int, output dimensions from each encoder
            n_classes: int, number of classification classes
            hidden_dim: int, hidden layer dimension
            dropout: float, dropout rate
        """
        super().__init__()
        total_dim = sum(encoder_dims)
        self.classifier = nn.Sequential(
            nn.Linear(total_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, n_classes),
        )

    def forward(self, embeddings):
        """
        Args:
            embeddings: list of tensors, each (batch, encoder_dim_i)
        Returns:
            logits: (batch, n_classes)
        """
        concat = torch.cat(embeddings, dim=-1)
        return self.classifier(concat)


# ===========================================================================
# Save importable architecture module
# ===========================================================================
def save_architecture_module():
    """Save the encoder/fusion classes as a standalone Python module."""
    module_path = os.path.join(OUTDIR, "multimodal_architecture.py")

    module_code = '''#!/usr/bin/env python3
"""
multimodal_architecture.py
Importable multi-modal encoder and fusion model definitions.

Generated by 72_multimodal_hooks.py.
Usage:
    from multimodal_architecture import (
        BulkRNAEncoder, AtlasFeatureEncoder, ScRNAEncoder,
        ATACEncoder, ProteomicsEncoder, FusionModel,
    )
"""

import torch
import torch.nn as nn


class MultiModalEncoder(nn.Module):
    """Base class for modality-specific encoders."""

    def __init__(self, output_dim: int):
        super().__init__()
        self.output_dim = output_dim

    def encode(self, x):
        raise NotImplementedError

    def forward(self, x):
        return self.encode(x)


class BulkRNAEncoder(MultiModalEncoder):
    """Uses pre-trained VAE encoder to map 3000 genes -> 64-dim embedding."""

    def __init__(self, n_genes=3000, latent_dim=64, dropout=0.3, freeze=True):
        super().__init__(output_dim=latent_dim)
        self.encoder = nn.Sequential(
            nn.Linear(n_genes, 2048),
            nn.BatchNorm1d(2048),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(2048, 512),
            nn.BatchNorm1d(512),
            nn.ReLU(),
            nn.Dropout(dropout),
        )
        self.fc_mu = nn.Linear(512, latent_dim)
        self.freeze = freeze

    def load_from_vae(self, vae_state_dict):
        """Load encoder weights from a trained ConditionalVAE state dict."""
        mapping = {}
        for key in vae_state_dict:
            if key.startswith("encoder."):
                mapping[key] = vae_state_dict[key]
            elif key in ("fc_mu.weight", "fc_mu.bias"):
                mapping[key] = vae_state_dict[key]
        self.load_state_dict(mapping, strict=False)
        if self.freeze:
            for param in self.parameters():
                param.requires_grad = False

    def encode(self, x):
        h = self.encoder(x)
        return self.fc_mu(h)


class AtlasFeatureEncoder(MultiModalEncoder):
    """Maps N atlas numeric columns -> 32-dim embedding."""

    def __init__(self, n_features, output_dim=32, dropout=0.2):
        super().__init__(output_dim=output_dim)
        self.net = nn.Sequential(
            nn.Linear(n_features, 64),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(64, output_dim),
        )

    def encode(self, x):
        return self.net(x)


class ScRNAEncoder(MultiModalEncoder):
    """Placeholder: aggregated cell-type proportion/expression embeddings."""

    def __init__(self, n_celltypes=10, output_dim=32):
        super().__init__(output_dim=output_dim)
        self.net = nn.Sequential(
            nn.Linear(n_celltypes, 64),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(64, output_dim),
        )

    def encode(self, x):
        return self.net(x)


class ATACEncoder(MultiModalEncoder):
    """Placeholder: chromatin accessibility peak features."""

    def __init__(self, n_peaks=500, output_dim=32):
        super().__init__(output_dim=output_dim)
        self.net = nn.Sequential(
            nn.Linear(n_peaks, 256),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(256, output_dim),
        )

    def encode(self, x):
        return self.net(x)


class ProteomicsEncoder(MultiModalEncoder):
    """Placeholder: protein abundance features."""

    def __init__(self, n_proteins=200, output_dim=32):
        super().__init__(output_dim=output_dim)
        self.net = nn.Sequential(
            nn.Linear(n_proteins, 128),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(128, output_dim),
        )

    def encode(self, x):
        return self.net(x)


class FusionModel(nn.Module):
    """
    Concatenates modality embeddings -> shared MLP -> prediction head.

    Args:
        encoder_dims: list of int, output dimensions from each encoder
        n_classes: int, number of classification classes
        hidden_dim: int, hidden layer dimension
        dropout: float, dropout rate
    """

    def __init__(self, encoder_dims, n_classes, hidden_dim=64, dropout=0.3):
        super().__init__()
        total_dim = sum(encoder_dims)
        self.classifier = nn.Sequential(
            nn.Linear(total_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, n_classes),
        )

    def forward(self, embeddings):
        """
        Args:
            embeddings: list of tensors, each (batch, encoder_dim_i)
        Returns:
            logits: (batch, n_classes)
        """
        concat = torch.cat(embeddings, dim=-1)
        return self.classifier(concat)
'''

    with open(module_path, "w") as f:
        f.write(module_code)
    print(f"  Architecture module saved: {module_path}")


# ===========================================================================
# Load data
# ===========================================================================
print("Loading data...")

with h5py.File(H5_PATH, "r") as h5:
    expr = h5["rank_expression"][:].T.astype(np.float32)  # samples x genes
    gene_names = [g.decode() for g in h5["gene_names"][:]]
    sample_ids = [s.decode() for s in h5["sample_ids"][:]]

    meta = {}
    for key in h5["metadata"].keys():
        raw = h5["metadata"][key][:]
        if raw.dtype.kind == "S":
            meta[key] = np.array([v.decode() for v in raw])
        else:
            meta[key] = raw.astype(np.float32)

    # Atlas features
    atlas_features = None
    atlas_genes = None
    atlas_feature_names = None
    if "atlas_features" in h5:
        atlas_features = h5["atlas_features/features"][:].astype(np.float32)
        atlas_genes = [g.decode() for g in h5["atlas_features/gene"][:]]
        atlas_feature_names = [n.decode() for n in h5["atlas_features/feature_names"][:]]
        print(f"  Atlas features: {atlas_features.shape}")
    else:
        print("  WARNING: No atlas features found in HDF5. Fusion demo will be limited.")

N_SAMPLES, N_GENES = expr.shape
fib_stage = meta["fib_stage"].astype(np.int64)
loco_folds = meta["loco_fold_fibrosis"]
unique_folds = sorted(set(f for f in loco_folds if f not in ("excluded", "NA")))

print(f"  Samples: {N_SAMPLES}, Genes: {N_GENES}")
print(f"  Fibrosis-labeled: {(fib_stage >= 0).sum()}")
print(f"  LOCO folds: {len(unique_folds)}")


# ===========================================================================
# Build per-sample atlas feature vectors
# ===========================================================================
def build_sample_atlas_features(expr_genes, atlas_features, atlas_genes, atlas_feature_names):
    """
    Build a per-sample atlas feature vector by matching atlas genes to expression genes.

    For each sample, the atlas features describe the *genes* in the expression panel.
    We aggregate by mean to create a fixed-size per-sample feature vector.
    Also select individual atlas columns that vary across genes as direct features.

    Returns: (n_samples, n_atlas_features) array
    """
    if atlas_features is None:
        return None, None

    # Build gene -> atlas row index
    atlas_gene_to_idx = {g: i for i, g in enumerate(atlas_genes)}

    # Match expression genes to atlas rows
    matched_idx = []
    for gene in expr_genes:
        if gene in atlas_gene_to_idx:
            matched_idx.append(atlas_gene_to_idx[gene])
        else:
            matched_idx.append(-1)

    matched_idx = np.array(matched_idx)
    n_matched = (matched_idx >= 0).sum()
    print(f"  Atlas gene matching: {n_matched}/{len(expr_genes)} expression genes found in atlas")

    if n_matched == 0:
        return None, None

    # Strategy: For each atlas column, compute a per-sample weighted sum
    # using the expression values as weights. This creates a
    # "expression-weighted atlas profile" per sample.
    # Filter to atlas columns with sufficient variance
    col_var = atlas_features.var(axis=0)
    informative_cols = np.where(col_var > 1e-6)[0]
    print(f"  Informative atlas columns: {len(informative_cols)} / {len(atlas_feature_names)}")

    if len(informative_cols) == 0:
        return None, None

    # Limit to top 50 most variable columns to keep feature dim manageable
    if len(informative_cols) > 50:
        top_var_idx = np.argsort(col_var[informative_cols])[-50:]
        informative_cols = informative_cols[top_var_idx]

    sel_atlas = atlas_features[:, informative_cols]
    sel_names = [atlas_feature_names[i] for i in informative_cols]

    # Build gene-level atlas matrix aligned with expression matrix
    gene_atlas = np.zeros((len(expr_genes), len(informative_cols)), dtype=np.float32)
    for i, aidx in enumerate(matched_idx):
        if aidx >= 0:
            gene_atlas[i] = sel_atlas[aidx]

    # Per-sample feature: expression-weighted mean of atlas features
    # expr is (n_samples, n_genes), gene_atlas is (n_genes, n_atlas_cols)
    # Normalize expression weights to sum to 1 per sample
    expr_norm = expr / (expr.sum(axis=1, keepdims=True) + 1e-8)
    sample_atlas = expr_norm @ gene_atlas  # (n_samples, n_atlas_cols)

    # Standardize
    mu = sample_atlas.mean(axis=0)
    sd = sample_atlas.std(axis=0) + 1e-8
    sample_atlas = (sample_atlas - mu) / sd

    print(f"  Per-sample atlas features: {sample_atlas.shape}")
    return sample_atlas, sel_names


sample_atlas, atlas_col_names = build_sample_atlas_features(
    gene_names, atlas_features, atlas_genes, atlas_feature_names,
)


# ===========================================================================
# Multi-modal Dataset
# ===========================================================================
class MultiModalDataset(Dataset):
    """Dataset providing expression + atlas features + fibrosis labels."""

    def __init__(self, expr, atlas, fib_stage, indices=None):
        self.indices = indices if indices is not None else np.arange(len(expr))
        self.expr = torch.tensor(expr[self.indices], dtype=torch.float32)
        self.fib = torch.tensor(fib_stage[self.indices], dtype=torch.long)
        if atlas is not None:
            self.atlas = torch.tensor(atlas[self.indices], dtype=torch.float32)
        else:
            self.atlas = None

    def __len__(self):
        return len(self.indices)

    def __getitem__(self, idx):
        if self.atlas is not None:
            return self.expr[idx], self.atlas[idx], self.fib[idx]
        return self.expr[idx], torch.zeros(1), self.fib[idx]


# ===========================================================================
# Compute class weights
# ===========================================================================
def compute_class_weights(y, n_classes):
    """Balanced class weights for cross-entropy."""
    valid = y >= 0
    if valid.sum() == 0:
        return torch.ones(n_classes)
    y_valid = y[valid]
    classes, counts = np.unique(y_valid, return_counts=True)
    weights = torch.ones(n_classes)
    for cls, cnt in zip(classes, counts):
        if int(cls) < n_classes:
            weights[int(cls)] = len(y_valid) / (n_classes * cnt)
    return weights


# ===========================================================================
# Train and evaluate a fusion or single-modality model
# ===========================================================================
def train_and_eval(bulk_encoder, atlas_encoder, fusion_model,
                   train_loader, val_loader, test_loader,
                   fib_weights, n_epochs=100, lr=1e-4, patience=15,
                   mode="fusion"):
    """
    Train fusion model (or single-modality) and return test metrics.

    mode: "fusion" uses both bulk+atlas, "bulk_only" uses bulk only,
          "atlas_only" uses atlas only.
    """
    # Collect trainable parameters
    params = []
    if bulk_encoder is not None and mode != "atlas_only":
        params.extend([p for p in bulk_encoder.parameters() if p.requires_grad])
    if atlas_encoder is not None and mode != "bulk_only":
        params.extend(atlas_encoder.parameters())
    params.extend(fusion_model.parameters())

    if len(params) == 0:
        print("    WARNING: No trainable parameters")
        return {}

    optimizer = torch.optim.Adam(params, lr=lr)
    fib_weights = fib_weights.to(DEVICE)

    best_val_loss = float("inf")
    best_states = {}
    patience_counter = 0

    for epoch in range(1, n_epochs + 1):
        # Train
        if bulk_encoder is not None:
            bulk_encoder.train()
        if atlas_encoder is not None:
            atlas_encoder.train()
        fusion_model.train()

        for x_batch, a_batch, f_batch in train_loader:
            x_batch = x_batch.to(DEVICE)
            a_batch = a_batch.to(DEVICE)
            f_batch = f_batch.to(DEVICE)

            embeddings = []
            if mode in ("fusion", "bulk_only") and bulk_encoder is not None:
                embeddings.append(bulk_encoder(x_batch))
            if mode in ("fusion", "atlas_only") and atlas_encoder is not None:
                embeddings.append(atlas_encoder(a_batch))

            if len(embeddings) == 0:
                continue

            logits = fusion_model(embeddings)

            fib_mask = f_batch >= 0
            if fib_mask.sum() == 0:
                continue

            loss = F.cross_entropy(logits[fib_mask], f_batch[fib_mask], weight=fib_weights)

            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(params, max_norm=5.0)
            optimizer.step()

        # Validation
        if bulk_encoder is not None:
            bulk_encoder.eval()
        if atlas_encoder is not None:
            atlas_encoder.eval()
        fusion_model.eval()

        val_loss = 0.0
        val_n = 0
        with torch.no_grad():
            for x_batch, a_batch, f_batch in val_loader:
                x_batch = x_batch.to(DEVICE)
                a_batch = a_batch.to(DEVICE)
                f_batch = f_batch.to(DEVICE)

                embeddings = []
                if mode in ("fusion", "bulk_only") and bulk_encoder is not None:
                    embeddings.append(bulk_encoder(x_batch))
                if mode in ("fusion", "atlas_only") and atlas_encoder is not None:
                    embeddings.append(atlas_encoder(a_batch))

                if len(embeddings) == 0:
                    continue

                logits = fusion_model(embeddings)
                fib_mask = f_batch >= 0
                if fib_mask.sum() > 0:
                    val_loss += F.cross_entropy(
                        logits[fib_mask], f_batch[fib_mask], weight=fib_weights,
                    ).item()
                    val_n += 1

        val_loss /= max(val_n, 1)

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_states = {}
            if bulk_encoder is not None:
                best_states["bulk"] = {k: v.cpu().clone()
                                        for k, v in bulk_encoder.state_dict().items()}
            if atlas_encoder is not None:
                best_states["atlas"] = {k: v.cpu().clone()
                                         for k, v in atlas_encoder.state_dict().items()}
            best_states["fusion"] = {k: v.cpu().clone()
                                      for k, v in fusion_model.state_dict().items()}
            patience_counter = 0
        else:
            patience_counter += 1

        if patience_counter >= patience:
            break

    # Restore best state
    if "bulk" in best_states and bulk_encoder is not None:
        bulk_encoder.load_state_dict(best_states["bulk"])
        bulk_encoder.to(DEVICE)
    if "atlas" in best_states and atlas_encoder is not None:
        atlas_encoder.load_state_dict(best_states["atlas"])
        atlas_encoder.to(DEVICE)
    if "fusion" in best_states:
        fusion_model.load_state_dict(best_states["fusion"])
        fusion_model.to(DEVICE)

    # Test evaluation
    if bulk_encoder is not None:
        bulk_encoder.eval()
    if atlas_encoder is not None:
        atlas_encoder.eval()
    fusion_model.eval()

    all_preds = []
    all_probs = []
    all_trues = []

    with torch.no_grad():
        for x_batch, a_batch, f_batch in test_loader:
            x_batch = x_batch.to(DEVICE)
            a_batch = a_batch.to(DEVICE)

            embeddings = []
            if mode in ("fusion", "bulk_only") and bulk_encoder is not None:
                embeddings.append(bulk_encoder(x_batch))
            if mode in ("fusion", "atlas_only") and atlas_encoder is not None:
                embeddings.append(atlas_encoder(a_batch))

            if len(embeddings) == 0:
                continue

            logits = fusion_model(embeddings)
            probs = F.softmax(logits, dim=1)
            preds = logits.argmax(dim=1)

            fib_mask = f_batch >= 0
            all_preds.extend(preds[fib_mask].cpu().numpy())
            all_probs.append(probs[fib_mask].cpu().numpy())
            all_trues.extend(f_batch[fib_mask].numpy())

    if len(all_trues) < 5:
        return {"balanced_accuracy": np.nan, "kappa": np.nan, "auroc": np.nan,
                "n_test": len(all_trues)}

    all_trues = np.array(all_trues)
    all_preds = np.array(all_preds)
    all_probs = np.vstack(all_probs)

    metrics = {
        "balanced_accuracy": balanced_accuracy_score(all_trues, all_preds),
        "kappa": cohen_kappa_score(all_trues, all_preds, weights="quadratic"),
        "n_test": len(all_trues),
    }

    try:
        metrics["auroc"] = roc_auc_score(
            all_trues, all_probs, multi_class="ovr", average="weighted",
        )
    except ValueError:
        metrics["auroc"] = np.nan

    return metrics


# ===========================================================================
# Main
# ===========================================================================
def main():
    t0 = time.time()

    # Save the importable architecture module
    save_architecture_module()

    # Check prerequisites
    if not os.path.exists(VAE_PATH):
        print(f"ERROR: VAE model not found at {VAE_PATH}")
        print("  Run 69_vae_training.py first.")
        sys.exit(1)

    has_atlas = sample_atlas is not None
    n_atlas_dim = sample_atlas.shape[1] if has_atlas else 0

    # Load VAE weights
    print("\nLoading pre-trained VAE encoder...")
    vae_state = torch.load(VAE_PATH, map_location="cpu", weights_only=True)
    print(f"  VAE state keys: {len(vae_state)}")

    # Fibrosis class weights
    fib_weights = compute_class_weights(fib_stage, 5)

    # -----------------------------------------------------------------------
    # LOCO Fusion Demo
    # -----------------------------------------------------------------------
    print(f"\n--- LOCO Fusion Demo ({len(unique_folds)} folds) ---")
    if has_atlas:
        print(f"  Atlas feature dim: {n_atlas_dim}")
    else:
        print("  No atlas features available; running bulk-only baseline")

    all_results = []

    for fold_name in unique_folds:
        print(f"\n  Fold: {fold_name}")

        fold_test_mask = np.array([f == fold_name for f in loco_folds])
        fold_train_all = ~fold_test_mask

        # Valid fibrosis samples only for train/val split
        valid = fib_stage >= 0
        train_cands = np.where(fold_train_all)[0]
        np.random.seed(42)
        perm = np.random.permutation(train_cands)
        n_tr = int(0.8 * len(perm))
        train_idx = perm[:n_tr]
        val_idx = perm[n_tr:]
        test_idx = np.where(fold_test_mask)[0]

        train_ds = MultiModalDataset(expr, sample_atlas, fib_stage, indices=train_idx)
        val_ds = MultiModalDataset(expr, sample_atlas, fib_stage, indices=val_idx)
        test_ds = MultiModalDataset(expr, sample_atlas, fib_stage, indices=test_idx)

        train_loader = DataLoader(train_ds, batch_size=32, shuffle=True,
                                  num_workers=2, pin_memory=True)
        val_loader = DataLoader(val_ds, batch_size=32, shuffle=False,
                                num_workers=2, pin_memory=True)
        test_loader = DataLoader(test_ds, batch_size=64, shuffle=False,
                                 num_workers=2, pin_memory=True)

        # ----- Experiment 1: Bulk RNA only -----
        bulk_enc = BulkRNAEncoder(n_genes=N_GENES, latent_dim=64, freeze=True).to(DEVICE)
        bulk_enc.load_from_vae(vae_state)

        fusion_bulk_only = FusionModel(
            encoder_dims=[64], n_classes=5, hidden_dim=64, dropout=0.3,
        ).to(DEVICE)

        m_bulk = train_and_eval(
            bulk_enc, None, fusion_bulk_only,
            train_loader, val_loader, test_loader,
            fib_weights, n_epochs=100, lr=1e-4, patience=15,
            mode="bulk_only",
        )
        m_bulk["fold"] = fold_name
        m_bulk["mode"] = "bulk_only"
        all_results.append(m_bulk)
        print(f"    Bulk-only:  bacc={m_bulk.get('balanced_accuracy', np.nan):.3f}  "
              f"kappa={m_bulk.get('kappa', np.nan):.3f}  "
              f"auroc={m_bulk.get('auroc', np.nan):.3f}")

        # ----- Experiment 2: Fusion (bulk + atlas) -----
        if has_atlas:
            bulk_enc_fuse = BulkRNAEncoder(
                n_genes=N_GENES, latent_dim=64, freeze=True,
            ).to(DEVICE)
            bulk_enc_fuse.load_from_vae(vae_state)

            atlas_enc = AtlasFeatureEncoder(
                n_features=n_atlas_dim, output_dim=32, dropout=0.2,
            ).to(DEVICE)

            fusion_model = FusionModel(
                encoder_dims=[64, 32], n_classes=5, hidden_dim=64, dropout=0.3,
            ).to(DEVICE)

            m_fuse = train_and_eval(
                bulk_enc_fuse, atlas_enc, fusion_model,
                train_loader, val_loader, test_loader,
                fib_weights, n_epochs=100, lr=1e-4, patience=15,
                mode="fusion",
            )
            m_fuse["fold"] = fold_name
            m_fuse["mode"] = "fusion"
            all_results.append(m_fuse)
            print(f"    Fusion:     bacc={m_fuse.get('balanced_accuracy', np.nan):.3f}  "
                  f"kappa={m_fuse.get('kappa', np.nan):.3f}  "
                  f"auroc={m_fuse.get('auroc', np.nan):.3f}")

            # ----- Experiment 3: Atlas only (ablation) -----
            atlas_enc_only = AtlasFeatureEncoder(
                n_features=n_atlas_dim, output_dim=32, dropout=0.2,
            ).to(DEVICE)

            fusion_atlas_only = FusionModel(
                encoder_dims=[32], n_classes=5, hidden_dim=64, dropout=0.3,
            ).to(DEVICE)

            m_atlas = train_and_eval(
                None, atlas_enc_only, fusion_atlas_only,
                train_loader, val_loader, test_loader,
                fib_weights, n_epochs=100, lr=1e-4, patience=15,
                mode="atlas_only",
            )
            m_atlas["fold"] = fold_name
            m_atlas["mode"] = "atlas_only"
            all_results.append(m_atlas)
            print(f"    Atlas-only: bacc={m_atlas.get('balanced_accuracy', np.nan):.3f}  "
                  f"kappa={m_atlas.get('kappa', np.nan):.3f}  "
                  f"auroc={m_atlas.get('auroc', np.nan):.3f}")

        # Free GPU
        torch.cuda.empty_cache()

    # -----------------------------------------------------------------------
    # Save results
    # -----------------------------------------------------------------------
    results_df = pd.DataFrame(all_results)
    results_path = os.path.join(OUTDIR, "fusion_demo_results.csv")
    results_df.to_csv(results_path, index=False)
    print(f"\nFusion results saved: {results_path}")

    # -----------------------------------------------------------------------
    # Modality contribution summary (ablation)
    # -----------------------------------------------------------------------
    print("\n--- Modality Contribution Summary ---")
    contribution_rows = []

    for metric in ["balanced_accuracy", "kappa", "auroc"]:
        bulk_vals = results_df[results_df["mode"] == "bulk_only"][metric].dropna().values
        if has_atlas:
            fuse_vals = results_df[results_df["mode"] == "fusion"][metric].dropna().values
            atlas_vals = results_df[results_df["mode"] == "atlas_only"][metric].dropna().values
        else:
            fuse_vals = np.array([])
            atlas_vals = np.array([])

        row = {
            "metric": metric,
            "bulk_only_mean": np.mean(bulk_vals) if len(bulk_vals) > 0 else np.nan,
            "bulk_only_std": np.std(bulk_vals) if len(bulk_vals) > 0 else np.nan,
        }
        if len(fuse_vals) > 0:
            row["fusion_mean"] = np.mean(fuse_vals)
            row["fusion_std"] = np.std(fuse_vals)
            row["fusion_delta"] = np.mean(fuse_vals) - np.mean(bulk_vals)
        if len(atlas_vals) > 0:
            row["atlas_only_mean"] = np.mean(atlas_vals)
            row["atlas_only_std"] = np.std(atlas_vals)

        contribution_rows.append(row)

        desc = f"  {metric}: bulk={row['bulk_only_mean']:.3f}"
        if "fusion_mean" in row:
            desc += f"  fusion={row['fusion_mean']:.3f} (delta={row['fusion_delta']:+.3f})"
        if "atlas_only_mean" in row:
            desc += f"  atlas={row['atlas_only_mean']:.3f}"
        print(desc)

    contrib_df = pd.DataFrame(contribution_rows)
    contrib_path = os.path.join(OUTDIR, "modality_contribution.csv")
    contrib_df.to_csv(contrib_path, index=False)
    print(f"\nModality contribution saved: {contrib_path}")

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    dt = time.time() - t0
    print(f"\n{'=' * 60}")
    print(f"72_multimodal_hooks.py completed in {dt / 60:.1f} min")
    print(f"  Architecture:  {os.path.join(OUTDIR, 'multimodal_architecture.py')}")
    print(f"  Fusion demo:   {results_path}")
    print(f"  Contribution:  {contrib_path}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
