#!/usr/bin/env python3
"""
90_nas_contrastive_vae.py
NAS-specific contrastive VAE — fork of Script 69 with NAS-focused ordinal
contrastive loss and multi-target contrastive weighting.

Architecture:
  Same ConditionalVAE as Script 69:
  Encoder: Linear(3000->2048)->BN->ReLU->Drop(0.3)->Linear(2048->512)->BN->ReLU->Drop(0.3)->(mu,logvar: 512->64)
  Decoder: Linear(64+8->512)->BN->ReLU->Drop(0.3)->Linear(512->2048)->BN->ReLU->Linear(2048->3000)
  Stage conditioning: Embedding(max_stages=10, dim=8) conditioned on NAS (not fibrosis)
  Domain discriminator: Linear(64->32)->ReLU->Linear(32->n_datasets) with gradient reversal

Key changes from Script 69:
  1. Contrastive loss uses NAS score (0-8) — positive: |NAS_diff|<=1, negative: |NAS_diff|>=2
  2. Multi-target contrastive: 0.7 NAS + 0.3 fibrosis (retain disease structure)
  3. Stage embedding conditioned on NAS (not fibrosis)
  4. Semi-supervised: 1,444 samples for recon+KL+domain; 660 NAS-annotated for NAS contrastive;
     1,128 fibrosis-annotated for fibrosis contrastive
  5. After VAE training, fit classical models on 64-dim NAS embeddings for NAS prediction

Training: Adam lr=1e-4, batch=64, 200 epochs, early stopping patience=30.
KL warmup: beta 0->1 over 50 epochs.
LOCO: 5-fold NAS LOCO + full model.

Input:  results/staging_classifier/prepared_data.h5
Output: results/staging_classifier/
  - nas_vae_model_full.pt
  - nas_embeddings_all_samples.csv (64-dim, 1444 samples)
  - nas_vae_training_curves.csv
  - nas_vae_classical_results.csv
  - nas_vae_model_summary.csv

SLURM: gpu partition, 1xL40S, 8 CPUs, 64GB RAM, 48h
Env:    micromamba activate rapids_singlecell
"""

import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
import h5py
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------
SEED = 42


def set_seed(seed=SEED):
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


set_seed()

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
INT = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
OUTDIR = os.path.join(INT, "results/staging_classifier")
LOCO_DIR = os.path.join(OUTDIR, "nas_embeddings_loco")
LOCO_GENE_DIR = os.path.join(OUTDIR, "loco_gene_lists")
os.makedirs(LOCO_DIR, exist_ok=True)

H5_PATH = os.path.join(OUTDIR, "prepared_data.h5")
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

print("=" * 70)
print("90: NAS-Specific Contrastive VAE with Multi-Target Loss")
print("=" * 70)
print(f"Device:    {DEVICE}")
print(f"Input:     {H5_PATH}")
print(f"Output:    {OUTDIR}")
print()

# ---------------------------------------------------------------------------
# Load HDF5 data
# ---------------------------------------------------------------------------
def load_h5_data(h5_path):
    """Load expression matrix and metadata from prepared_data.h5."""
    with h5py.File(h5_path, "r") as h5:
        # Expression: genes x samples -> transpose to samples x genes
        expr = h5["rank_expression"][:].T.astype(np.float32)
        gene_names = [g.decode() for g in h5["gene_names"][:]]
        sample_ids = [s.decode() for s in h5["sample_ids"][:]]

        # Metadata
        meta = {}
        for key in h5["metadata"].keys():
            raw = h5["metadata"][key][:]
            if raw.dtype.kind == "S":
                meta[key] = np.array([v.decode() for v in raw])
            else:
                meta[key] = raw.astype(np.float32)

    return {
        "expr": expr,
        "gene_names": gene_names,
        "sample_ids": sample_ids,
        "meta": meta,
    }


print("Loading data...")
data = load_h5_data(H5_PATH)
N_SAMPLES, N_GENES = data["expr"].shape
print(f"  Samples: {N_SAMPLES}, Genes: {N_GENES}")

# Build gene name -> index mapping for fold-specific subsetting
_gene_name_to_idx_90 = {g: i for i, g in enumerate(data["gene_names"])}

# Load global feature candidates as fallback
_global_fc_path_90 = os.path.join(OUTDIR, "feature_candidates_3000.csv")
if os.path.exists(_global_fc_path_90):
    _global_fc_90 = pd.read_csv(_global_fc_path_90)
    _feature_genes_global_90 = _global_fc_90["gene"].tolist()
else:
    _feature_genes_global_90 = list(data["gene_names"])


def load_fold_gene_indices_90(fold_id):
    """Return column indices into expr for the per-fold gene list.

    Falls back to the global feature_candidates_3000.csv if the fold-specific
    file does not exist.
    """
    fold_gene_file = os.path.join(
        LOCO_GENE_DIR, f"fold_{fold_id}_feature_candidates_3000.csv")
    if os.path.exists(fold_gene_file):
        fold_genes = pd.read_csv(fold_gene_file)["gene"].tolist()
        print(f"    Using per-fold gene list for fold: {fold_id}")
    else:
        fold_genes = _feature_genes_global_90
        print(f"    Fold gene list not found for {fold_id} -- using global fallback")
    return np.array([_gene_name_to_idx_90[g] for g in fold_genes
                     if g in _gene_name_to_idx_90], dtype=np.int64)

# Build dataset-to-index mapping
datasets = data["meta"]["dataset"]
unique_datasets = sorted(set(datasets))
dataset_to_idx = {d: i for i, d in enumerate(unique_datasets)}
dataset_idx = np.array([dataset_to_idx[d] for d in datasets], dtype=np.int64)
N_DATASETS = len(unique_datasets)
print(f"  Datasets: {N_DATASETS} ({', '.join(unique_datasets)})")

# Build NAS score array (-1 = unlabeled)
nas_score = data["meta"]["nas_score"].astype(np.int64)
n_nas_labeled = (nas_score >= 0).sum()
print(f"  NAS-labeled:      {n_nas_labeled} / {N_SAMPLES}")

# Build fibrosis stage array (-1 = unlabeled)
fib_stage = data["meta"]["fib_stage"].astype(np.int64)
n_fib_labeled = (fib_stage >= 0).sum()
print(f"  Fibrosis-labeled: {n_fib_labeled} / {N_SAMPLES}")

# Build NAS 4-group array (-1 = unlabeled)
nas_group4 = data["meta"]["nas_group4"].astype(np.int64)

# Build NAS 3-class: 0 (NAS 0-2), 1 (NAS 3-5), 2 (NAS 6-8)
nas_3class = np.full(N_SAMPLES, -1, dtype=np.int64)
nas_valid_mask = nas_score >= 0
nas_3class[nas_valid_mask & (nas_score <= 2)] = 0
nas_3class[nas_valid_mask & (nas_score >= 3) & (nas_score <= 5)] = 1
nas_3class[nas_valid_mask & (nas_score >= 6)] = 2
print(f"  NAS 3-class distribution: {dict(zip(*np.unique(nas_3class[nas_3class >= 0], return_counts=True)))}")

# Build NAS>=5 binary
nas_ge5 = data["meta"]["nas_ge5"].astype(np.int64)

# NAS distribution
nas_vals, nas_counts = np.unique(nas_score[nas_score >= 0], return_counts=True)
for v, c in zip(nas_vals, nas_counts):
    print(f"    NAS {v}: n={c}")

# LOCO fold assignments (NAS)
loco_folds_nas = data["meta"]["loco_fold_nas"]
unique_folds_nas = sorted(set(f for f in loco_folds_nas if f != "excluded"))
print(f"  LOCO NAS folds: {len(unique_folds_nas)} ({', '.join(unique_folds_nas)})")


# ---------------------------------------------------------------------------
# PyTorch Dataset (multi-target)
# ---------------------------------------------------------------------------
class NASDataset(Dataset):
    """Expression dataset with NAS, fibrosis, and dataset labels."""

    def __init__(self, expr, nas_score, fib_stage, dataset_idx, indices=None):
        self.indices = indices if indices is not None else np.arange(len(expr))
        self.expr = torch.tensor(expr[self.indices], dtype=torch.float32)
        self.nas_score = torch.tensor(nas_score[self.indices], dtype=torch.long)
        self.fib_stage = torch.tensor(fib_stage[self.indices], dtype=torch.long)
        self.dataset_idx = torch.tensor(dataset_idx[self.indices], dtype=torch.long)

    def __len__(self):
        return len(self.indices)

    def __getitem__(self, idx):
        return (self.expr[idx], self.nas_score[idx],
                self.fib_stage[idx], self.dataset_idx[idx])


# ---------------------------------------------------------------------------
# Gradient Reversal Layer
# ---------------------------------------------------------------------------
class GradientReversalFunction(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, lmbda):
        ctx.lmbda = lmbda
        return x.clone()

    @staticmethod
    def backward(ctx, grad_output):
        return -ctx.lmbda * grad_output, None


class GradientReversalLayer(nn.Module):
    def __init__(self, lmbda=0.1):
        super().__init__()
        self.lmbda = lmbda

    def forward(self, x):
        return GradientReversalFunction.apply(x, self.lmbda)


# ---------------------------------------------------------------------------
# Conditional VAE (NAS-conditioned)
# ---------------------------------------------------------------------------
class NASConditionalVAE(nn.Module):
    """
    Conditional VAE with NAS-focused ordinal contrastive + domain adversarial.
    Stage embedding conditioned on NAS score (not fibrosis).
    """

    def __init__(self, n_genes=3000, latent_dim=64, stage_embed_dim=8,
                 max_stages=10, n_datasets=10, dropout=0.3):
        super().__init__()
        self.latent_dim = latent_dim
        self.stage_embed_dim = stage_embed_dim

        # Encoder
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
        self.fc_logvar = nn.Linear(512, latent_dim)

        # NAS stage embedding for conditioning the decoder (NAS 0-9)
        self.stage_embed = nn.Embedding(max_stages, stage_embed_dim)

        # Decoder
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim + stage_embed_dim, 512),
            nn.BatchNorm1d(512),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(512, 2048),
            nn.BatchNorm1d(2048),
            nn.ReLU(),
            nn.Linear(2048, n_genes),
        )

        # Domain discriminator with gradient reversal
        self.grl = GradientReversalLayer(lmbda=0.1)
        self.domain_head = nn.Sequential(
            nn.Linear(latent_dim, 32),
            nn.ReLU(),
            nn.Linear(32, n_datasets),
        )

    def encode(self, x):
        h = self.encoder(x)
        return self.fc_mu(h), self.fc_logvar(h)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z, nas_labels):
        # For unlabeled samples (NAS=-1), use default embedding (index 0)
        safe_labels = nas_labels.clone()
        safe_labels[safe_labels < 0] = 0
        s_embed = self.stage_embed(safe_labels)
        zs = torch.cat([z, s_embed], dim=-1)
        return self.decoder(zs)

    def forward(self, x, nas_labels):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        x_recon = self.decode(z, nas_labels)
        domain_logits = self.domain_head(self.grl(mu))  # deterministic for domain
        return x_recon, mu, logvar, z, domain_logits


# ---------------------------------------------------------------------------
# Ordinal Contrastive Loss (triplet with ordinal distance weighting)
# ---------------------------------------------------------------------------
def ordinal_contrastive_loss(z, stages, base_margin=0.5):
    """
    Triplet margin loss with ordinal distance weighting.
    Anchor: sample at stage k
    Positive: sample at stage k +/- 1
    Negative: sample at stage >= k+2 away
    Margin weighted by |stage_anchor - stage_negative|

    Only uses samples with valid staging (stage >= 0).
    Returns loss (scalar) or None if no valid triplets.
    """
    mask = stages >= 0
    if mask.sum() < 3:
        return None

    z_labeled = z[mask]
    s_labeled = stages[mask]
    n = z_labeled.shape[0]

    anchor_idx, pos_idx, neg_idx, margins = [], [], [], []
    stage_vals = s_labeled.cpu().numpy()

    for a in range(n):
        sa = stage_vals[a]
        pos_cands = np.where(np.abs(stage_vals - sa) <= 1)[0]
        neg_cands = np.where(np.abs(stage_vals - sa) >= 2)[0]

        # Exclude self from positives
        pos_cands = pos_cands[pos_cands != a]

        if len(pos_cands) == 0 or len(neg_cands) == 0:
            continue

        # Sample one positive and one negative per anchor
        p = pos_cands[np.random.randint(len(pos_cands))]
        ne = neg_cands[np.random.randint(len(neg_cands))]
        m = base_margin * abs(sa - stage_vals[ne])

        anchor_idx.append(a)
        pos_idx.append(p)
        neg_idx.append(ne)
        margins.append(m)

    if len(anchor_idx) == 0:
        return None

    anchor_idx = torch.tensor(anchor_idx, device=z_labeled.device)
    pos_idx = torch.tensor(pos_idx, device=z_labeled.device)
    neg_idx = torch.tensor(neg_idx, device=z_labeled.device)
    margins_t = torch.tensor(margins, dtype=torch.float32, device=z_labeled.device)

    z_a = z_labeled[anchor_idx]
    z_p = z_labeled[pos_idx]
    z_n = z_labeled[neg_idx]

    d_pos = F.pairwise_distance(z_a, z_p)
    d_neg = F.pairwise_distance(z_a, z_n)
    losses = F.relu(d_pos - d_neg + margins_t)
    return losses.mean()


# ---------------------------------------------------------------------------
# Multi-target contrastive: 0.7 NAS + 0.3 fibrosis
# ---------------------------------------------------------------------------
def multi_target_contrastive_loss(z, nas_scores, fib_stages,
                                  w_nas=0.7, w_fib=0.3, base_margin=0.5):
    """
    Weighted combination of NAS and fibrosis ordinal contrastive losses.
    NAS: uses nas_scores (0-8), only NAS-labeled samples
    Fibrosis: uses fib_stages (0-4), only fibrosis-labeled samples
    """
    loss_nas = ordinal_contrastive_loss(z, nas_scores, base_margin=base_margin)
    loss_fib = ordinal_contrastive_loss(z, fib_stages, base_margin=base_margin)

    total = torch.tensor(0.0, device=z.device)
    has_loss = False

    if loss_nas is not None:
        total = total + w_nas * loss_nas
        has_loss = True
    if loss_fib is not None:
        total = total + w_fib * loss_fib
        has_loss = True

    if not has_loss:
        return None, None, None

    return total, loss_nas, loss_fib


# ---------------------------------------------------------------------------
# Training loop
# ---------------------------------------------------------------------------
def train_nas_vae(model, train_loader, val_loader, n_epochs=200, lr=1e-4,
                  patience=30, kl_warmup=50, w_contrastive=0.5, w_domain=0.1,
                  w_nas=0.7, w_fib=0.3):
    """Train the NAS conditional VAE. Returns training curves and best model state."""

    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    best_val_recon = float("inf")
    best_state = None
    patience_counter = 0
    curves = []

    for epoch in range(1, n_epochs + 1):
        model.train()
        epoch_losses = {
            "recon": 0, "kl": 0, "contr_nas": 0, "contr_fib": 0,
            "contr_total": 0, "domain": 0, "total": 0,
        }
        n_batches = 0

        # KL annealing
        beta = min(1.0, epoch / kl_warmup) if kl_warmup > 0 else 1.0

        for x_batch, nas_batch, fib_batch, ds_batch in train_loader:
            x_batch = x_batch.to(DEVICE)
            nas_batch = nas_batch.to(DEVICE)
            fib_batch = fib_batch.to(DEVICE)
            ds_batch = ds_batch.to(DEVICE)

            x_recon, mu, logvar, z, domain_logits = model(x_batch, nas_batch)

            # 1. Reconstruction loss (MSE)
            loss_recon = F.mse_loss(x_recon, x_batch, reduction="mean")

            # 2. KL divergence
            loss_kl = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())

            # 3. Multi-target ordinal contrastive (NAS + fibrosis)
            loss_contr, loss_nas_c, loss_fib_c = multi_target_contrastive_loss(
                z, nas_batch, fib_batch, w_nas=w_nas, w_fib=w_fib,
            )
            if loss_contr is None:
                contr_val = 0.0
                nas_c_val = 0.0
                fib_c_val = 0.0
                loss_contr = torch.tensor(0.0, device=DEVICE)
            else:
                contr_val = loss_contr.item()
                nas_c_val = loss_nas_c.item() if loss_nas_c is not None else 0.0
                fib_c_val = loss_fib_c.item() if loss_fib_c is not None else 0.0

            # 4. Domain adversarial
            loss_domain = F.cross_entropy(domain_logits, ds_batch)

            # Total loss
            total = (loss_recon
                     + beta * loss_kl
                     + w_contrastive * loss_contr
                     + w_domain * loss_domain)

            optimizer.zero_grad()
            total.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()

            epoch_losses["recon"] += loss_recon.item()
            epoch_losses["kl"] += loss_kl.item()
            epoch_losses["contr_nas"] += nas_c_val
            epoch_losses["contr_fib"] += fib_c_val
            epoch_losses["contr_total"] += contr_val
            epoch_losses["domain"] += loss_domain.item()
            epoch_losses["total"] += total.item()
            n_batches += 1

        # Average over batches
        for k in epoch_losses:
            epoch_losses[k] /= max(n_batches, 1)

        # Validation (reconstruction loss only)
        model.eval()
        val_recon = 0.0
        val_batches = 0
        with torch.no_grad():
            for x_batch, nas_batch, fib_batch, ds_batch in val_loader:
                x_batch = x_batch.to(DEVICE)
                nas_batch = nas_batch.to(DEVICE)
                x_recon, mu, logvar, z, domain_logits = model(x_batch, nas_batch)
                val_recon += F.mse_loss(x_recon, x_batch, reduction="mean").item()
                val_batches += 1
        val_recon /= max(val_batches, 1)

        curves.append({
            "epoch": epoch,
            "beta": beta,
            "train_recon": epoch_losses["recon"],
            "train_kl": epoch_losses["kl"],
            "train_contr_nas": epoch_losses["contr_nas"],
            "train_contr_fib": epoch_losses["contr_fib"],
            "train_contr_total": epoch_losses["contr_total"],
            "train_domain": epoch_losses["domain"],
            "train_total": epoch_losses["total"],
            "val_recon": val_recon,
        })

        # Early stopping on validation reconstruction
        if val_recon < best_val_recon:
            best_val_recon = val_recon
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            patience_counter = 0
        else:
            patience_counter += 1

        if epoch % 20 == 0 or epoch == 1:
            print(f"  Epoch {epoch:3d} | recon={epoch_losses['recon']:.4f} "
                  f"kl={epoch_losses['kl']:.4f} "
                  f"c_nas={epoch_losses['contr_nas']:.4f} "
                  f"c_fib={epoch_losses['contr_fib']:.4f} "
                  f"dom={epoch_losses['domain']:.4f} | "
                  f"val_recon={val_recon:.4f} | "
                  f"beta={beta:.3f} | pat={patience_counter}/{patience}")

        if patience_counter >= patience:
            print(f"  Early stopping at epoch {epoch} (patience={patience})")
            break

    return best_state, pd.DataFrame(curves)


# ---------------------------------------------------------------------------
# Embed all samples with a trained model
# ---------------------------------------------------------------------------
@torch.no_grad()
def embed_samples(model, expr, nas_score, fib_stage, dataset_idx, batch_size=256):
    """Generate latent embeddings (mu) for all samples."""
    model.eval()
    ds = NASDataset(expr, nas_score, fib_stage, dataset_idx)
    loader = DataLoader(ds, batch_size=batch_size, shuffle=False)
    all_mu = []
    for x_batch, nas_batch, fib_batch, ds_batch in loader:
        x_batch = x_batch.to(DEVICE)
        nas_batch = nas_batch.to(DEVICE)
        mu, _ = model.encode(x_batch)
        all_mu.append(mu.cpu().numpy())
    return np.vstack(all_mu)


# ---------------------------------------------------------------------------
# Classical model evaluation on NAS embeddings
# ---------------------------------------------------------------------------
def evaluate_classical_models(embeddings_by_fold, fold_test_masks,
                              nas_score, fib_stage, nas_group4, nas_3class, nas_ge5,
                              loco_folds_nas, unique_folds):
    """
    Fit classical models on 64-dim NAS embeddings for NAS prediction tasks.
    Returns per-model per-fold results DataFrame and best model summary.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.svm import SVC
    from sklearn.metrics import (
        accuracy_score, roc_auc_score, f1_score, cohen_kappa_score,
        balanced_accuracy_score,
    )
    from sklearn.preprocessing import StandardScaler

    # Try to import xgboost
    try:
        from xgboost import XGBClassifier
        has_xgb = True
    except ImportError:
        print("  WARNING: xgboost not available, skipping XGB models")
        has_xgb = False

    # Define targets
    targets = {
        "nas_ge5_binary": {
            "labels": nas_ge5,
            "desc": "NAS>=5 binary",
            "multiclass": False,
        },
        "nas_group4": {
            "labels": nas_group4,
            "desc": "NAS 4-group ordinal",
            "multiclass": True,
        },
        "nas_3class": {
            "labels": nas_3class,
            "desc": "NAS 3-class (Low/Med/High)",
            "multiclass": True,
        },
        "nas_9class": {
            "labels": nas_score,
            "desc": "NAS 9-class ordinal (0-8)",
            "multiclass": True,
        },
    }

    # Define models
    def get_models(n_classes, multiclass):
        models = {
            "LogisticRegression": LogisticRegression(
                penalty="elasticnet", solver="saga", l1_ratio=0.5,
                max_iter=5000, random_state=SEED, C=1.0,
            ),
            "RandomForest": RandomForestClassifier(
                n_estimators=500, max_depth=10, min_samples_leaf=5,
                random_state=SEED, n_jobs=-1,
            ),
            "SVM": SVC(
                kernel="rbf", probability=True, random_state=SEED, C=1.0,
            ),
        }
        if has_xgb:
            models["XGBoost"] = XGBClassifier(
                n_estimators=500, max_depth=6, learning_rate=0.05,
                subsample=0.8, colsample_bytree=0.8,
                random_state=SEED, n_jobs=-1, verbosity=0,
                objective="multi:softprob" if multiclass else "binary:logistic",
                num_class=n_classes if multiclass else None,
                eval_metric="mlogloss" if multiclass else "logloss",
            )
        return models

    results = []

    for target_name, target_info in targets.items():
        labels = target_info["labels"]
        multiclass = target_info["multiclass"]
        valid_mask = labels >= 0
        n_classes = len(np.unique(labels[valid_mask]))

        print(f"\n  Target: {target_info['desc']} ({n_classes} classes, "
              f"{valid_mask.sum()} labeled samples)")

        models = get_models(n_classes, multiclass)

        for fold_name in unique_folds:
            fold_idx = list(unique_folds).index(fold_name)

            # Get embeddings for this fold's model
            emb = embeddings_by_fold[fold_name]

            # Test set: held-out fold samples with valid labels
            test_mask = fold_test_masks[fold_name] & valid_mask
            # Train set: all other NAS-labeled samples
            train_mask = ~fold_test_masks[fold_name] & valid_mask

            if test_mask.sum() < 5 or train_mask.sum() < 10:
                print(f"    Fold {fold_name}: skipped (train={train_mask.sum()}, "
                      f"test={test_mask.sum()})")
                continue

            X_train = emb[train_mask]
            y_train = labels[train_mask]
            X_test = emb[test_mask]
            y_test = labels[test_mask]

            # Standardize embeddings
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)

            for model_name, clf in models.items():
                try:
                    clf.fit(X_train, y_train)
                    y_pred = clf.predict(X_test)
                    y_prob = clf.predict_proba(X_test)

                    acc = accuracy_score(y_test, y_pred)
                    bal_acc = balanced_accuracy_score(y_test, y_pred)
                    f1 = f1_score(y_test, y_pred, average="macro", zero_division=0)
                    kappa = cohen_kappa_score(y_test, y_pred)

                    # AUROC
                    if multiclass:
                        try:
                            auroc = roc_auc_score(
                                y_test, y_prob, multi_class="ovr", average="macro",
                            )
                        except ValueError:
                            auroc = np.nan
                    else:
                        try:
                            auroc = roc_auc_score(y_test, y_prob[:, 1])
                        except (ValueError, IndexError):
                            auroc = np.nan

                    # Quadratic weighted kappa (for ordinal targets)
                    qwk = cohen_kappa_score(y_test, y_pred, weights="quadratic")

                    results.append({
                        "target": target_name,
                        "target_desc": target_info["desc"],
                        "fold": fold_name,
                        "model": model_name,
                        "n_train": int(train_mask.sum()),
                        "n_test": int(test_mask.sum()),
                        "n_classes": n_classes,
                        "accuracy": acc,
                        "balanced_accuracy": bal_acc,
                        "auroc": auroc,
                        "f1_macro": f1,
                        "kappa": kappa,
                        "qwk": qwk,
                    })

                    if fold_idx == 0:
                        print(f"    {model_name:25s} | acc={acc:.3f} "
                              f"bal_acc={bal_acc:.3f} auroc={auroc:.3f} "
                              f"f1={f1:.3f} qwk={qwk:.3f}")

                except Exception as e:
                    print(f"    {model_name} FAILED on fold {fold_name}: {e}")
                    results.append({
                        "target": target_name,
                        "target_desc": target_info["desc"],
                        "fold": fold_name,
                        "model": model_name,
                        "n_train": int(train_mask.sum()),
                        "n_test": int(test_mask.sum()),
                        "n_classes": n_classes,
                        "accuracy": np.nan,
                        "balanced_accuracy": np.nan,
                        "auroc": np.nan,
                        "f1_macro": np.nan,
                        "kappa": np.nan,
                        "qwk": np.nan,
                    })

    results_df = pd.DataFrame(results)

    # Compute summary: mean metrics per target x model
    if len(results_df) > 0:
        summary = (
            results_df
            .groupby(["target", "target_desc", "model"])
            .agg(
                mean_accuracy=("accuracy", "mean"),
                std_accuracy=("accuracy", "std"),
                mean_balanced_accuracy=("balanced_accuracy", "mean"),
                mean_auroc=("auroc", "mean"),
                std_auroc=("auroc", "std"),
                mean_f1_macro=("f1_macro", "mean"),
                mean_qwk=("qwk", "mean"),
                n_folds=("fold", "count"),
            )
            .reset_index()
            .sort_values(["target", "mean_auroc"], ascending=[True, False])
        )
    else:
        summary = pd.DataFrame()

    return results_df, summary


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()

    # -----------------------------------------------------------------------
    # Full model training (all samples)
    # -----------------------------------------------------------------------
    print("\n--- Training full NAS VAE (all samples) ---")

    # 80/20 train/val split (stratified by dataset)
    set_seed()
    train_mask = np.zeros(N_SAMPLES, dtype=bool)
    for ds_name in unique_datasets:
        ds_mask = np.where(datasets == ds_name)[0]
        n_train = max(int(0.8 * len(ds_mask)), 1)
        perm = np.random.permutation(ds_mask)
        train_mask[perm[:n_train]] = True
    val_mask = ~train_mask

    print(f"  Train: {train_mask.sum()}, Val: {val_mask.sum()}")

    train_ds = NASDataset(data["expr"], nas_score, fib_stage, dataset_idx,
                          indices=np.where(train_mask)[0])
    val_ds = NASDataset(data["expr"], nas_score, fib_stage, dataset_idx,
                        indices=np.where(val_mask)[0])
    train_loader = DataLoader(train_ds, batch_size=64, shuffle=True,
                              num_workers=2, pin_memory=True, drop_last=False)
    val_loader = DataLoader(val_ds, batch_size=64, shuffle=False,
                            num_workers=2, pin_memory=True)

    model_full = NASConditionalVAE(
        n_genes=N_GENES, latent_dim=64, stage_embed_dim=8,
        max_stages=10, n_datasets=N_DATASETS, dropout=0.3,
    ).to(DEVICE)

    n_params = sum(p.numel() for p in model_full.parameters())
    print(f"  Model params: {n_params:,}")

    best_state_full, curves_full = train_nas_vae(
        model_full, train_loader, val_loader,
        n_epochs=200, lr=1e-4, patience=30, kl_warmup=50,
        w_contrastive=0.5, w_domain=0.1, w_nas=0.7, w_fib=0.3,
    )

    # Save training curves
    curves_path = os.path.join(OUTDIR, "nas_vae_training_curves.csv")
    curves_full.to_csv(curves_path, index=False)
    print(f"  Training curves saved ({len(curves_full)} epochs)")

    # Load best state and embed all samples
    model_full.load_state_dict(best_state_full)
    model_full.to(DEVICE)
    embeddings_full = embed_samples(
        model_full, data["expr"], nas_score, fib_stage, dataset_idx,
    )

    emb_df = pd.DataFrame(
        embeddings_full,
        index=data["sample_ids"],
        columns=[f"z{i}" for i in range(64)],
    )
    emb_df.index.name = "sample_id"
    emb_path = os.path.join(OUTDIR, "nas_embeddings_all_samples.csv")
    emb_df.to_csv(emb_path)
    print(f"  Full embeddings: {emb_df.shape}")

    # Save model weights
    model_path = os.path.join(OUTDIR, "nas_vae_model_full.pt")
    torch.save(best_state_full, model_path)
    print(f"  Full model saved: {model_path}")

    # -----------------------------------------------------------------------
    # LOCO VAE training (per NAS fold) — 5 folds
    # -----------------------------------------------------------------------
    print(f"\n--- LOCO NAS VAE training ({len(unique_folds_nas)} folds) ---")

    embeddings_by_fold = {}
    fold_test_masks = {}

    for fold_name in unique_folds_nas:
        print(f"\n  Fold: {fold_name} (held out)")
        t_fold = time.time()

        # Load per-fold gene list and subset expression (removes leakage)
        fold_gene_idx = load_fold_gene_indices_90(fold_name)
        expr_fold = data["expr"][:, fold_gene_idx]
        n_genes_fold = expr_fold.shape[1]

        # Held-out = this fold; train = everything else
        fold_mask_test = np.array([f == fold_name for f in loco_folds_nas])
        fold_mask_train_all = ~fold_mask_test

        fold_test_masks[fold_name] = fold_mask_test

        # Further split train into train/val (80/20 within training set)
        train_candidates = np.where(fold_mask_train_all)[0]
        set_seed()
        perm = np.random.permutation(train_candidates)
        n_train = int(0.8 * len(perm))
        fold_train_idx = perm[:n_train]
        fold_val_idx = perm[n_train:]

        print(f"    Train: {len(fold_train_idx)}, Val: {len(fold_val_idx)}, "
              f"Test (held-out): {fold_mask_test.sum()}, Features: {n_genes_fold}")

        fold_train_ds = NASDataset(expr_fold, nas_score, fib_stage, dataset_idx,
                                   indices=fold_train_idx)
        fold_val_ds = NASDataset(expr_fold, nas_score, fib_stage, dataset_idx,
                                 indices=fold_val_idx)
        fold_train_loader = DataLoader(fold_train_ds, batch_size=64, shuffle=True,
                                       num_workers=2, pin_memory=True)
        fold_val_loader = DataLoader(fold_val_ds, batch_size=64, shuffle=False,
                                     num_workers=2, pin_memory=True)

        model_fold = NASConditionalVAE(
            n_genes=n_genes_fold, latent_dim=64, stage_embed_dim=8,
            max_stages=10, n_datasets=N_DATASETS, dropout=0.3,
        ).to(DEVICE)

        best_state_fold, _ = train_nas_vae(
            model_fold, fold_train_loader, fold_val_loader,
            n_epochs=200, lr=1e-4, patience=30, kl_warmup=50,
            w_contrastive=0.5, w_domain=0.1, w_nas=0.7, w_fib=0.3,
        )

        # Embed ALL samples with this fold's model (using fold-specific genes)
        model_fold.load_state_dict(best_state_fold)
        model_fold.to(DEVICE)
        emb_fold = embed_samples(
            model_fold, expr_fold, nas_score, fib_stage, dataset_idx,
        )

        embeddings_by_fold[fold_name] = emb_fold

        emb_fold_df = pd.DataFrame(
            emb_fold,
            index=data["sample_ids"],
            columns=[f"z{i}" for i in range(64)],
        )
        emb_fold_df.index.name = "sample_id"
        emb_fold_df["is_heldout"] = fold_mask_test.astype(int)
        fold_csv = os.path.join(LOCO_DIR, f"nas_embeddings_loco_{fold_name}.csv")
        emb_fold_df.to_csv(fold_csv)
        dt = time.time() - t_fold
        print(f"    Saved {fold_csv} ({dt:.0f}s)")

        # Free GPU memory between folds
        del model_fold
        torch.cuda.empty_cache()

    # -----------------------------------------------------------------------
    # Classical model evaluation on NAS embeddings
    # -----------------------------------------------------------------------
    print("\n--- Classical model evaluation on NAS embeddings ---")

    results_df, summary_df = evaluate_classical_models(
        embeddings_by_fold=embeddings_by_fold,
        fold_test_masks=fold_test_masks,
        nas_score=nas_score,
        fib_stage=fib_stage,
        nas_group4=nas_group4,
        nas_3class=nas_3class,
        nas_ge5=nas_ge5,
        loco_folds_nas=loco_folds_nas,
        unique_folds=unique_folds_nas,
    )

    # Save results
    results_path = os.path.join(OUTDIR, "nas_vae_classical_results.csv")
    results_df.to_csv(results_path, index=False)
    print(f"\n  Classical results saved: {results_path} ({len(results_df)} rows)")

    summary_path = os.path.join(OUTDIR, "nas_vae_model_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"  Model summary saved:    {summary_path} ({len(summary_df)} rows)")

    # Print summary table
    if len(summary_df) > 0:
        print("\n  ===== Best model per target =====")
        for target in summary_df["target"].unique():
            sub = summary_df[summary_df["target"] == target].iloc[0]
            print(f"    {sub['target_desc']:35s} | "
                  f"best={sub['model']:20s} | "
                  f"AUROC={sub['mean_auroc']:.3f}+/-{sub['std_auroc']:.3f} | "
                  f"bal_acc={sub['mean_balanced_accuracy']:.3f} | "
                  f"QWK={sub['mean_qwk']:.3f}")

    # -----------------------------------------------------------------------
    # UMAP visualization
    # -----------------------------------------------------------------------
    print("\n--- Generating UMAP visualization ---")
    try:
        from sklearn.decomposition import PCA
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Try GPU UMAP first, fall back to CPU
        try:
            from cuml.manifold import UMAP as cuUMAP
            reducer = cuUMAP(n_components=2, n_neighbors=30, min_dist=0.3,
                             random_state=42)
        except ImportError:
            reducer = None

        if reducer is not None:
            umap_coords = reducer.fit_transform(embeddings_full)
        else:
            # Fallback: PCA -> t-SNE
            pca = PCA(n_components=20, random_state=42)
            pca_coords = pca.fit_transform(embeddings_full)
            from sklearn.manifold import TSNE
            tsne = TSNE(n_components=2, perplexity=30, random_state=42)
            umap_coords = tsne.fit_transform(pca_coords)

        fig, axes = plt.subplots(1, 3, figsize=(21, 6))

        # Panel 1: Color by NAS score
        ax = axes[0]
        nas_labeled = nas_score >= 0
        nas_unlabeled = nas_score < 0
        if nas_unlabeled.any():
            ax.scatter(umap_coords[nas_unlabeled, 0], umap_coords[nas_unlabeled, 1],
                       c="lightgrey", s=8, alpha=0.3, label="Unlabeled")
        if nas_labeled.any():
            sc = ax.scatter(umap_coords[nas_labeled, 0], umap_coords[nas_labeled, 1],
                            c=nas_score[nas_labeled], cmap="RdYlBu_r", s=15, alpha=0.7,
                            vmin=0, vmax=8)
            plt.colorbar(sc, ax=ax, label="NAS Score")
        ax.set_title("Latent Space — NAS Score")
        ax.set_xlabel("Dim 1")
        ax.set_ylabel("Dim 2")

        # Panel 2: Color by fibrosis stage
        ax = axes[1]
        fib_labeled = fib_stage >= 0
        fib_unlabeled = fib_stage < 0
        if fib_unlabeled.any():
            ax.scatter(umap_coords[fib_unlabeled, 0], umap_coords[fib_unlabeled, 1],
                       c="lightgrey", s=8, alpha=0.3, label="Unlabeled")
        if fib_labeled.any():
            sc = ax.scatter(umap_coords[fib_labeled, 0], umap_coords[fib_labeled, 1],
                            c=fib_stage[fib_labeled], cmap="viridis", s=15, alpha=0.7,
                            vmin=0, vmax=4)
            plt.colorbar(sc, ax=ax, label="Fibrosis Stage")
        ax.set_title("Latent Space — Fibrosis Stage")
        ax.set_xlabel("Dim 1")
        ax.set_ylabel("Dim 2")

        # Panel 3: Color by dataset (batch)
        ax = axes[2]
        for i, ds_name in enumerate(unique_datasets):
            mask = datasets == ds_name
            ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                       s=10, alpha=0.5, label=ds_name)
        ax.set_title("Latent Space — Dataset (Batch)")
        ax.set_xlabel("Dim 1")
        ax.set_ylabel("Dim 2")
        ax.legend(fontsize=5, ncol=2, loc="best")

        plt.tight_layout()
        fig_path = os.path.join(OUTDIR, "nas_latent_space_umap.png")
        plt.savefig(fig_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  UMAP saved: {fig_path}")
    except Exception as e:
        print(f"  UMAP generation failed (non-fatal): {e}")

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    dt_total = time.time() - t0
    print(f"\n{'=' * 70}")
    print(f"90_nas_contrastive_vae.py completed in {dt_total / 60:.1f} min")
    print(f"  Full embeddings: {emb_path}")
    print(f"  LOCO embeddings: {LOCO_DIR}/ ({len(unique_folds_nas)} folds)")
    print(f"  Model:           {model_path}")
    print(f"  Training curves: {curves_path}")
    print(f"  Classical:       {results_path}")
    print(f"  Summary:         {summary_path}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
