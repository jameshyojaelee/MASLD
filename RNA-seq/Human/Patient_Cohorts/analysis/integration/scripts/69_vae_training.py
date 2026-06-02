#!/usr/bin/env python3
"""
69_vae_training.py
Conditional VAE with ordinal contrastive loss and domain adaptation.

Architecture:
  Encoder: Linear(3000->2048)->BN->ReLU->Drop(0.3)->Linear(2048->512)->BN->ReLU->Drop(0.3)->(mu,logvar: 512->64)
  Decoder: Linear(64+8->512)->BN->ReLU->Drop(0.3)->Linear(512->2048)->BN->ReLU->Linear(2048->3000)
  Stage conditioning: Embedding(max_stages=10, dim=8) concatenated with z
  Domain discriminator: Linear(64->32)->ReLU->Linear(32->n_datasets) with gradient reversal

Losses:
  1. Reconstruction (MSE, weight=1.0)
  2. KL divergence (beta-VAE warmup 0->1 over 50 epochs)
  3. Ordinal contrastive (triplet margin, weight=0.5) -- staged samples only
  4. Domain adversarial (CE, weight=0.1, gradient reversal lambda=0.1)

Semi-supervised: unlabeled samples contribute to recon+KL only.

Training: Adam lr=1e-4, batch=64, 200 epochs, early stopping patience=30.
LOCO: per-fibrosis-fold VAE + full VAE for final embeddings.

Input:  results/staging_classifier/prepared_data.h5
Output: results/staging_classifier/
  - vae_model_full.pt
  - embeddings_all_samples.csv (64-dim, sample_id index)
  - embeddings_loco/*.csv
  - vae_training_curves.csv
  - latent_space_umap.png

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
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
INT = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
OUTDIR = os.path.join(INT, "results/staging_classifier")
LOCO_DIR = os.path.join(OUTDIR, "embeddings_loco")
LOCO_GENE_DIR = os.path.join(OUTDIR, "loco_gene_lists")
os.makedirs(LOCO_DIR, exist_ok=True)

H5_PATH = os.path.join(OUTDIR, "prepared_data.h5")
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

print("=" * 60)
print("69: Conditional VAE with Ordinal Contrastive + Domain Adapt")
print("=" * 60)
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

        # Atlas features (for Script 72)
        atlas_features = None
        atlas_genes = None
        atlas_feature_names = None
        if "atlas_features" in h5:
            atlas_features = h5["atlas_features/features"][:].astype(np.float32)
            atlas_genes = [g.decode() for g in h5["atlas_features/gene"][:]]
            atlas_feature_names = [n.decode() for n in h5["atlas_features/feature_names"][:]]

    return {
        "expr": expr,
        "gene_names": gene_names,
        "sample_ids": sample_ids,
        "meta": meta,
        "atlas_features": atlas_features,
        "atlas_genes": atlas_genes,
        "atlas_feature_names": atlas_feature_names,
    }


print("Loading data...")
data = load_h5_data(H5_PATH)
N_SAMPLES, N_GENES = data["expr"].shape
print(f"  Samples: {N_SAMPLES}, Genes: {N_GENES}")

# Build gene name -> index mapping for fold-specific subsetting
_gene_name_to_idx = {g: i for i, g in enumerate(data["gene_names"])}

# Load global feature candidates as fallback
_global_fc_path = os.path.join(OUTDIR, "feature_candidates_3000.csv")
if os.path.exists(_global_fc_path):
    _global_fc = pd.read_csv(_global_fc_path)
    _feature_genes_global = _global_fc["gene"].tolist()
else:
    _feature_genes_global = list(data["gene_names"])


def load_fold_gene_indices_69(fold_id):
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
        fold_genes = _feature_genes_global
        print(f"    Fold gene list not found for {fold_id} -- using global fallback")
    return np.array([_gene_name_to_idx[g] for g in fold_genes
                     if g in _gene_name_to_idx], dtype=np.int64)

# Build dataset-to-index mapping
datasets = data["meta"]["dataset"]
unique_datasets = sorted(set(datasets))
dataset_to_idx = {d: i for i, d in enumerate(unique_datasets)}
dataset_idx = np.array([dataset_to_idx[d] for d in datasets], dtype=np.int64)
N_DATASETS = len(unique_datasets)
print(f"  Datasets: {N_DATASETS} ({', '.join(unique_datasets)})")

# Build fibrosis stage array (-1 = unlabeled)
fib_stage = data["meta"]["fib_stage"].astype(np.int64)
n_labeled = (fib_stage >= 0).sum()
print(f"  Fibrosis-labeled: {n_labeled} / {N_SAMPLES}")

# LOCO fold assignments
loco_folds = data["meta"]["loco_fold_fibrosis"]
unique_folds = sorted(set(f for f in loco_folds if f != "excluded"))
print(f"  LOCO fibrosis folds: {len(unique_folds)} ({', '.join(unique_folds)})")


# ---------------------------------------------------------------------------
# PyTorch Dataset
# ---------------------------------------------------------------------------
class MASLDDataset(Dataset):
    """Expression dataset with staging and dataset labels."""

    def __init__(self, expr, fib_stage, dataset_idx, indices=None):
        self.indices = indices if indices is not None else np.arange(len(expr))
        self.expr = torch.tensor(expr[self.indices], dtype=torch.float32)
        self.fib_stage = torch.tensor(fib_stage[self.indices], dtype=torch.long)
        self.dataset_idx = torch.tensor(dataset_idx[self.indices], dtype=torch.long)

    def __len__(self):
        return len(self.indices)

    def __getitem__(self, idx):
        return self.expr[idx], self.fib_stage[idx], self.dataset_idx[idx]


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
# Conditional VAE
# ---------------------------------------------------------------------------
class ConditionalVAE(nn.Module):
    """Conditional VAE with ordinal contrastive + domain adversarial heads."""

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

        # Stage embedding for conditioning the decoder
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

    def decode(self, z, stage_labels):
        # For unlabeled samples (stage=-1), use a default embedding (index 0)
        safe_labels = stage_labels.clone()
        safe_labels[safe_labels < 0] = 0
        s_embed = self.stage_embed(safe_labels)
        zs = torch.cat([z, s_embed], dim=-1)
        return self.decoder(zs)

    def forward(self, x, stage_labels):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        x_recon = self.decode(z, stage_labels)
        domain_logits = self.domain_head(self.grl(mu))  # Use mu (deterministic) for domain
        return x_recon, mu, logvar, z, domain_logits


# ---------------------------------------------------------------------------
# Ordinal Contrastive Loss (triplet with ordinal distance weighting)
# ---------------------------------------------------------------------------
def ordinal_contrastive_loss(z, stages, base_margin=0.5):
    """
    Triplet margin loss with ordinal distance weighting.
    Only uses samples with valid staging (stage >= 0).
    Returns loss (scalar) or None if no valid triplets in this batch.
    """
    # Filter to labeled samples
    mask = stages >= 0
    if mask.sum() < 3:
        return None

    z_labeled = z[mask]
    s_labeled = stages[mask]
    n = z_labeled.shape[0]

    # Build all valid triplets: anchor, positive (|da-dp| <= 1), negative (|da-dn| >= 2)
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

        # Sample one positive and one negative per anchor (efficient)
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

    # Triplet margin loss: max(0, d(a,p) - d(a,n) + margin)
    d_pos = F.pairwise_distance(z_a, z_p)
    d_neg = F.pairwise_distance(z_a, z_n)
    losses = F.relu(d_pos - d_neg + margins_t)
    return losses.mean()


# ---------------------------------------------------------------------------
# Training loop
# ---------------------------------------------------------------------------
def train_vae(model, train_loader, val_loader, n_epochs=200, lr=1e-4,
              patience=30, kl_warmup=50, w_contrastive=0.5, w_domain=0.1):
    """Train the conditional VAE. Returns training curves and best model state."""

    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    best_val_recon = float("inf")
    best_state = None
    patience_counter = 0
    curves = []

    for epoch in range(1, n_epochs + 1):
        model.train()
        epoch_losses = {"recon": 0, "kl": 0, "contrastive": 0, "domain": 0, "total": 0}
        n_batches = 0

        # KL annealing: beta goes from 0 to 1 over first kl_warmup epochs
        beta = min(1.0, epoch / kl_warmup) if kl_warmup > 0 else 1.0

        for x_batch, stage_batch, ds_batch in train_loader:
            x_batch = x_batch.to(DEVICE)
            stage_batch = stage_batch.to(DEVICE)
            ds_batch = ds_batch.to(DEVICE)

            x_recon, mu, logvar, z, domain_logits = model(x_batch, stage_batch)

            # 1. Reconstruction loss (MSE)
            loss_recon = F.mse_loss(x_recon, x_batch, reduction="mean")

            # 2. KL divergence
            loss_kl = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())

            # 3. Ordinal contrastive (staged samples only)
            loss_contrastive = ordinal_contrastive_loss(z, stage_batch)
            if loss_contrastive is None:
                loss_contrastive_val = 0.0
                loss_contrastive = torch.tensor(0.0, device=DEVICE)
            else:
                loss_contrastive_val = loss_contrastive.item()

            # 4. Domain adversarial
            loss_domain = F.cross_entropy(domain_logits, ds_batch)

            # Total loss
            total = (loss_recon
                     + beta * loss_kl
                     + w_contrastive * loss_contrastive
                     + w_domain * loss_domain)

            optimizer.zero_grad()
            total.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()

            epoch_losses["recon"] += loss_recon.item()
            epoch_losses["kl"] += loss_kl.item()
            epoch_losses["contrastive"] += loss_contrastive_val
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
            for x_batch, stage_batch, ds_batch in val_loader:
                x_batch = x_batch.to(DEVICE)
                stage_batch = stage_batch.to(DEVICE)
                x_recon, mu, logvar, z, domain_logits = model(x_batch, stage_batch)
                val_recon += F.mse_loss(x_recon, x_batch, reduction="mean").item()
                val_batches += 1
        val_recon /= max(val_batches, 1)

        curves.append({
            "epoch": epoch,
            "beta": beta,
            "train_recon": epoch_losses["recon"],
            "train_kl": epoch_losses["kl"],
            "train_contrastive": epoch_losses["contrastive"],
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
                  f"kl={epoch_losses['kl']:.4f} contr={epoch_losses['contrastive']:.4f} "
                  f"dom={epoch_losses['domain']:.4f} | val_recon={val_recon:.4f} "
                  f"| beta={beta:.3f} | pat={patience_counter}/{patience}")

        if patience_counter >= patience:
            print(f"  Early stopping at epoch {epoch} (patience={patience})")
            break

    return best_state, pd.DataFrame(curves)


# ---------------------------------------------------------------------------
# Embed all samples with a trained model
# ---------------------------------------------------------------------------
@torch.no_grad()
def embed_samples(model, expr, fib_stage, dataset_idx, batch_size=256):
    """Generate latent embeddings (mu) for all samples."""
    model.eval()
    ds = MASLDDataset(expr, fib_stage, dataset_idx)
    loader = DataLoader(ds, batch_size=batch_size, shuffle=False)
    all_mu = []
    for x_batch, stage_batch, ds_batch in loader:
        x_batch = x_batch.to(DEVICE)
        stage_batch = stage_batch.to(DEVICE)
        mu, _ = model.encode(x_batch)
        all_mu.append(mu.cpu().numpy())
    return np.vstack(all_mu)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()

    # -----------------------------------------------------------------------
    # Full model training (all samples)
    # -----------------------------------------------------------------------
    print("\n--- Training full VAE (all samples) ---")

    # 80/20 train/val split (stratified by dataset)
    np.random.seed(42)
    train_mask = np.zeros(N_SAMPLES, dtype=bool)
    for ds_name in unique_datasets:
        ds_mask = np.where(datasets == ds_name)[0]
        n_train = max(int(0.8 * len(ds_mask)), 1)
        perm = np.random.permutation(ds_mask)
        train_mask[perm[:n_train]] = True
    val_mask = ~train_mask

    print(f"  Train: {train_mask.sum()}, Val: {val_mask.sum()}")

    train_ds = MASLDDataset(data["expr"], fib_stage, dataset_idx,
                            indices=np.where(train_mask)[0])
    val_ds = MASLDDataset(data["expr"], fib_stage, dataset_idx,
                          indices=np.where(val_mask)[0])
    train_loader = DataLoader(train_ds, batch_size=64, shuffle=True,
                              num_workers=2, pin_memory=True, drop_last=False)
    val_loader = DataLoader(val_ds, batch_size=64, shuffle=False,
                            num_workers=2, pin_memory=True)

    model_full = ConditionalVAE(
        n_genes=N_GENES, latent_dim=64, stage_embed_dim=8,
        max_stages=10, n_datasets=N_DATASETS, dropout=0.3,
    ).to(DEVICE)

    print(f"  Model params: {sum(p.numel() for p in model_full.parameters()):,}")

    best_state_full, curves_full = train_vae(
        model_full, train_loader, val_loader,
        n_epochs=200, lr=1e-4, patience=30, kl_warmup=50,
    )

    # Save training curves
    curves_full.to_csv(os.path.join(OUTDIR, "vae_training_curves.csv"), index=False)
    print(f"  Training curves saved ({len(curves_full)} epochs)")

    # Load best state and embed all samples
    model_full.load_state_dict(best_state_full)
    model_full.to(DEVICE)
    embeddings_full = embed_samples(model_full, data["expr"], fib_stage, dataset_idx)

    emb_df = pd.DataFrame(
        embeddings_full,
        index=data["sample_ids"],
        columns=[f"z{i}" for i in range(64)],
    )
    emb_df.index.name = "sample_id"
    emb_df.to_csv(os.path.join(OUTDIR, "embeddings_all_samples.csv"))
    print(f"  Full embeddings: {emb_df.shape}")

    # Save model weights
    torch.save(best_state_full, os.path.join(OUTDIR, "vae_model_full.pt"))
    print("  Full model saved")

    # -----------------------------------------------------------------------
    # LOCO VAE training (per fibrosis fold)
    # -----------------------------------------------------------------------
    print(f"\n--- LOCO VAE training ({len(unique_folds)} folds) ---")

    for fold_name in unique_folds:
        print(f"\n  Fold: {fold_name} (held out)")
        t_fold = time.time()

        # Load per-fold gene list and subset expression (removes leakage)
        fold_gene_idx = load_fold_gene_indices_69(fold_name)
        expr_fold = data["expr"][:, fold_gene_idx]
        n_genes_fold = expr_fold.shape[1]

        # Held-out = this fold; train = everything else
        fold_mask_test = np.array([f == fold_name for f in loco_folds])
        fold_mask_train_all = ~fold_mask_test

        # Further split train into train/val (80/20 within training set)
        train_candidates = np.where(fold_mask_train_all)[0]
        np.random.seed(42)
        perm = np.random.permutation(train_candidates)
        n_train = int(0.8 * len(perm))
        fold_train_idx = perm[:n_train]
        fold_val_idx = perm[n_train:]

        fold_train_ds = MASLDDataset(expr_fold, fib_stage, dataset_idx,
                                     indices=fold_train_idx)
        fold_val_ds = MASLDDataset(expr_fold, fib_stage, dataset_idx,
                                   indices=fold_val_idx)
        fold_train_loader = DataLoader(fold_train_ds, batch_size=64, shuffle=True,
                                       num_workers=2, pin_memory=True)
        fold_val_loader = DataLoader(fold_val_ds, batch_size=64, shuffle=False,
                                     num_workers=2, pin_memory=True)

        model_fold = ConditionalVAE(
            n_genes=n_genes_fold, latent_dim=64, stage_embed_dim=8,
            max_stages=10, n_datasets=N_DATASETS, dropout=0.3,
        ).to(DEVICE)

        best_state_fold, _ = train_vae(
            model_fold, fold_train_loader, fold_val_loader,
            n_epochs=200, lr=1e-4, patience=30, kl_warmup=50,
        )

        # Embed ALL samples with this fold's model (using fold-specific genes)
        model_fold.load_state_dict(best_state_fold)
        model_fold.to(DEVICE)
        emb_fold = embed_samples(model_fold, expr_fold, fib_stage, dataset_idx)

        emb_fold_df = pd.DataFrame(
            emb_fold,
            index=data["sample_ids"],
            columns=[f"z{i}" for i in range(64)],
        )
        emb_fold_df.index.name = "sample_id"
        emb_fold_df["is_heldout"] = fold_mask_test.astype(int)
        fold_csv = os.path.join(LOCO_DIR, f"embeddings_loco_{fold_name}.csv")
        emb_fold_df.to_csv(fold_csv)
        dt = time.time() - t_fold
        print(f"    Saved {fold_csv} ({dt:.0f}s)")

        # Free GPU memory between folds
        del model_fold
        torch.cuda.empty_cache()

    # -----------------------------------------------------------------------
    # UMAP visualization
    # -----------------------------------------------------------------------
    print("\n--- Generating UMAP ---")
    try:
        from sklearn.decomposition import PCA
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Try GPU UMAP first, fall back to sklearn
        try:
            from cuml.manifold import UMAP as cuUMAP
            reducer = cuUMAP(n_components=2, n_neighbors=30, min_dist=0.3,
                             random_state=42)
        except ImportError:
            from sklearn.manifold import TSNE
            # Use PCA to 2D as a simpler fallback
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

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Panel 1: Color by fibrosis stage
        ax = axes[0]
        labeled = fib_stage >= 0
        unlabeled = fib_stage < 0

        if unlabeled.any():
            ax.scatter(umap_coords[unlabeled, 0], umap_coords[unlabeled, 1],
                       c="lightgrey", s=8, alpha=0.3, label="Unlabeled")
        if labeled.any():
            sc = ax.scatter(umap_coords[labeled, 0], umap_coords[labeled, 1],
                            c=fib_stage[labeled], cmap="viridis", s=15, alpha=0.7,
                            vmin=0, vmax=4)
            plt.colorbar(sc, ax=ax, label="Fibrosis Stage")
        ax.set_title("Latent Space — Fibrosis Stage")
        ax.set_xlabel("Dim 1")
        ax.set_ylabel("Dim 2")

        # Panel 2: Color by dataset
        ax = axes[1]
        for i, ds_name in enumerate(unique_datasets):
            mask = datasets == ds_name
            ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                       s=10, alpha=0.5, label=ds_name)
        ax.set_title("Latent Space — Dataset (Batch)")
        ax.set_xlabel("Dim 1")
        ax.set_ylabel("Dim 2")
        ax.legend(fontsize=6, ncol=2, loc="best")

        plt.tight_layout()
        fig_path = os.path.join(OUTDIR, "latent_space_umap.png")
        plt.savefig(fig_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  UMAP saved: {fig_path}")
    except Exception as e:
        print(f"  UMAP generation failed (non-fatal): {e}")

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    dt_total = time.time() - t0
    print(f"\n{'=' * 60}")
    print(f"69_vae_training.py completed in {dt_total / 60:.1f} min")
    print(f"  Full embeddings: {os.path.join(OUTDIR, 'embeddings_all_samples.csv')}")
    print(f"  LOCO embeddings: {LOCO_DIR}/ ({len(unique_folds)} folds)")
    print(f"  Model:           {os.path.join(OUTDIR, 'vae_model_full.pt')}")
    print(f"  Training curves: {os.path.join(OUTDIR, 'vae_training_curves.csv')}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
