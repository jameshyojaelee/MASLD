#!/usr/bin/env python3
"""
93_fibrosis_transfer_nas.py
Test whether fibrosis-pretrained VAE features transfer to NAS prediction.

Two approaches:
  (A) Linear probe: Load fibrosis VAE embeddings (64-dim), train NAS ordinal
      models (elastic net, RF, XGBoost) on top. 5-fold NAS LOCO.
  (B) Fine-tune: Load fibrosis VAE weights, replace contrastive loss with
      NAS ordinal contrastive, fine-tune entire model for 50 epochs with
      low lr (1e-5 encoder, 1e-3 new head). 5-fold NAS LOCO.

Targets: NAS>=5 AUROC, NAS 4-group QWK, NAS 3-class accuracy.

Input:
  - results/staging_classifier/vae_model_full.pt (fibrosis-pretrained VAE)
  - results/staging_classifier/embeddings_all_samples.csv (fibrosis VAE embeddings)
  - results/staging_classifier/prepared_data.h5 (expression + NAS labels)

Output:
  - results/staging_classifier/fibrosis_transfer_results.csv
  - results/staging_classifier/transfer_vs_scratch_comparison.csv

SLURM: gpu partition, 1xL40S, 8 CPUs, 64GB RAM, 48h
Env:    micromamba activate rapids_singlecell
"""

import os
import sys
import time
import copy
import warnings
import numpy as np
import pandas as pd
import h5py
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    roc_auc_score, balanced_accuracy_score, cohen_kappa_score,
    f1_score, mean_absolute_error, accuracy_score,
)
import xgboost as xgb

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
EMB_PATH = os.path.join(OUTDIR, "embeddings_all_samples.csv")
VAE_PATH = os.path.join(OUTDIR, "vae_model_full.pt")

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

print("=" * 60)
print("93: Fibrosis-Pretrained VAE Transfer to NAS Prediction")
print("=" * 60)
print(f"Device:  {DEVICE}")
print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print()


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
print("Loading data...")

# Fibrosis VAE embeddings
emb_df = pd.read_csv(EMB_PATH, index_col="sample_id")
emb_cols = [c for c in emb_df.columns if c.startswith("z")]
print(f"  Fibrosis VAE embeddings: {emb_df.shape[0]} samples x {len(emb_cols)} dims")

# Expression + metadata from HDF5
with h5py.File(H5_PATH, "r") as h5:
    raw_expr = h5["rank_expression"][:].T.astype(np.float32)  # samples x genes
    gene_names = [g.decode() for g in h5["gene_names"][:]]
    sample_ids = [s.decode() for s in h5["sample_ids"][:]]

    meta = {}
    for key in h5["metadata"].keys():
        raw = h5["metadata"][key][:]
        if raw.dtype.kind == "S":
            meta[key] = np.array([v.decode() for v in raw])
        else:
            meta[key] = raw.astype(np.float32)

N_SAMPLES, N_GENES = raw_expr.shape
print(f"  Expression: {N_SAMPLES} samples x {N_GENES} genes")

# Build metadata
meta_df = pd.DataFrame(meta, index=sample_ids)
meta_df.index.name = "sample_id"

# Modeling metadata (for NAS 3-class)
modeling_meta = pd.read_csv(
    os.path.join(OUTDIR, "modeling_metadata.csv"), index_col="sample_id"
)

# Alignment check
assert list(emb_df.index) == sample_ids, "Embedding/expression sample order mismatch"

# NAS labels
nas_group = meta_df["nas_group"].values.astype(int)        # raw NAS score, -1 = missing
nas_group4 = meta_df["nas_group4"].values.astype(int)      # 4-group, -1 = missing
nas_ge5 = meta_df["nas_ge5"].values.astype(int)            # binary NAS>=5, -1 = missing
nas_3class = modeling_meta.loc[sample_ids, "nas_3class"].values.astype(int)  # 3-class

# Fibrosis labels (for fine-tuning reference)
fib_stage = meta_df["fib_stage"].values.astype(int)

# Dataset labels
datasets = meta_df["dataset"].values
unique_datasets = sorted(set(datasets))
dataset_to_idx = {d: i for i, d in enumerate(unique_datasets)}
dataset_idx = np.array([dataset_to_idx[d] for d in datasets], dtype=np.int64)
N_DATASETS = len(unique_datasets)

# LOCO NAS folds
loco_nas = meta_df["loco_fold_nas"].values
unique_nas_folds = sorted(set(f for f in loco_nas if f != "excluded" and f != "NA"))
print(f"  NAS LOCO folds: {len(unique_nas_folds)} ({', '.join(unique_nas_folds)})")

n_nas_labeled = (nas_group >= 0).sum()
print(f"  NAS-labeled: {n_nas_labeled} / {N_SAMPLES}")


# ---------------------------------------------------------------------------
# Target definitions
# ---------------------------------------------------------------------------
TARGETS = {
    "nas_ge5": {
        "y": nas_ge5,
        "n_classes": 2,
        "task": "binary",
        "description": "NAS>=5 (definite NASH)",
        "metric": "auroc",
    },
    "nas_4group": {
        "y": nas_group4,
        "n_classes": 4,
        "task": "multiclass",
        "description": "NAS 4-group",
        "metric": "qwk",
    },
    "nas_3class": {
        "y": nas_3class,
        "n_classes": 3,
        "task": "multiclass",
        "description": "NAS 3-class (Low/Med/High)",
        "metric": "accuracy",
    },
}


# ---------------------------------------------------------------------------
# Metrics helpers
# ---------------------------------------------------------------------------
def safe_auroc(y_true, y_prob, n_classes):
    """Compute AUROC handling edge cases."""
    try:
        if n_classes == 2:
            if y_prob.ndim == 2:
                return roc_auc_score(y_true, y_prob[:, 1])
            return roc_auc_score(y_true, y_prob)
        else:
            return roc_auc_score(y_true, y_prob, multi_class="ovr", average="weighted")
    except (ValueError, IndexError):
        return np.nan


def compute_metrics(y_true, y_pred, y_prob, n_classes):
    """Compute all relevant metrics."""
    metrics = {
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "accuracy": accuracy_score(y_true, y_pred),
        "f1_macro": f1_score(y_true, y_pred, average="macro", zero_division=0),
    }
    if n_classes > 2:
        metrics["qwk"] = cohen_kappa_score(y_true, y_pred, weights="quadratic")
        metrics["mae"] = mean_absolute_error(y_true, y_pred)
    else:
        metrics["qwk"] = cohen_kappa_score(y_true, y_pred)
        metrics["mae"] = np.nan
    if y_prob is not None:
        metrics["auroc"] = safe_auroc(y_true, y_prob, n_classes)
    else:
        metrics["auroc"] = np.nan
    return metrics


def compute_sample_weights(y):
    """Balanced sample weights for XGBoost."""
    classes, counts = np.unique(y, return_counts=True)
    n = len(y)
    nc = len(classes)
    weights = np.zeros(n)
    for cls, cnt in zip(classes, counts):
        weights[y == cls] = n / (nc * cnt)
    return weights


# =========================================================================
# APPROACH A: Linear probe on frozen fibrosis-VAE embeddings
# =========================================================================
def run_linear_probe():
    """Train NAS classifiers on frozen fibrosis-VAE embeddings, NAS LOCO CV."""
    print("\n" + "=" * 60)
    print("APPROACH A: Linear Probe on Fibrosis-VAE Embeddings")
    print("=" * 60)

    X_emb = emb_df[emb_cols].values.astype(np.float32)
    results = []

    for target_name, tinfo in TARGETS.items():
        y = tinfo["y"]
        n_classes = tinfo["n_classes"]
        valid = y >= 0

        print(f"\n  Target: {target_name} ({tinfo['description']}, "
              f"{valid.sum()} labeled samples)")

        for fold_name in unique_nas_folds:
            test_mask = (loco_nas == fold_name) & valid
            train_mask = (loco_nas != fold_name) & (loco_nas != "excluded") & \
                         (loco_nas != "NA") & valid

            if test_mask.sum() < 5 or train_mask.sum() < 10:
                continue

            y_train = y[train_mask].astype(int)
            y_test = y[test_mask].astype(int)

            if len(np.unique(y_train)) < 2:
                continue

            X_train = X_emb[train_mask]
            X_test = X_emb[test_mask]

            # Elastic net logistic
            models = {
                "ElasticNet": LogisticRegression(
                    penalty="elasticnet", l1_ratio=0.5, solver="saga",
                    class_weight="balanced", max_iter=5000, random_state=42,
                    C=1.0, n_jobs=-1,
                ),
                "RandomForest": RandomForestClassifier(
                    n_estimators=500, class_weight="balanced",
                    random_state=42, n_jobs=-1,
                ),
                "XGBoost": xgb.XGBClassifier(
                    n_estimators=300, learning_rate=0.05, max_depth=6,
                    random_state=42, n_jobs=-1,
                    eval_metric="mlogloss" if n_classes > 2 else "logloss",
                ),
            }

            for model_name, model in models.items():
                scaler = StandardScaler()
                X_tr_sc = scaler.fit_transform(X_train)
                X_te_sc = scaler.transform(X_test)

                fit_kwargs = {}
                if model_name == "XGBoost":
                    fit_kwargs["sample_weight"] = compute_sample_weights(y_train)

                try:
                    model.fit(X_tr_sc, y_train, **fit_kwargs)
                    y_pred = model.predict(X_te_sc)
                    try:
                        y_prob = model.predict_proba(X_te_sc)
                    except Exception:
                        y_prob = None

                    m = compute_metrics(y_test, y_pred, y_prob, n_classes)
                    m.update({
                        "approach": "A_linear_probe",
                        "target": target_name,
                        "model": model_name,
                        "fold": fold_name,
                        "n_train": int(train_mask.sum()),
                        "n_test": int(test_mask.sum()),
                    })
                    results.append(m)
                except Exception as e:
                    print(f"    {model_name} fold={fold_name}: {e}")

        # Print per-target summary
        target_res = [r for r in results if r["target"] == target_name
                      and r["approach"] == "A_linear_probe"]
        if target_res:
            df_t = pd.DataFrame(target_res)
            for mn in df_t["model"].unique():
                sub = df_t[df_t["model"] == mn]
                primary_metric = tinfo["metric"]
                val = sub[primary_metric].mean()
                print(f"    {mn}: mean {primary_metric}={val:.3f} "
                      f"({len(sub)} folds)")

    return results


# =========================================================================
# APPROACH B: Fine-tune fibrosis VAE for NAS
# =========================================================================

# Redefine the VAE architecture (must match Script 69)
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


class ConditionalVAE(nn.Module):
    """Conditional VAE matching Script 69 architecture."""

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
        safe_labels = stage_labels.clone()
        safe_labels[safe_labels < 0] = 0
        s_embed = self.stage_embed(safe_labels)
        zs = torch.cat([z, s_embed], dim=-1)
        return self.decoder(zs)

    def forward(self, x, stage_labels):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        x_recon = self.decode(z, stage_labels)
        domain_logits = self.domain_head(self.grl(mu))
        return x_recon, mu, logvar, z, domain_logits


# NAS prediction head (attached to VAE encoder)
class NASHead(nn.Module):
    """Prediction head for NAS targets."""

    def __init__(self, latent_dim=64, n_classes=4):
        super().__init__()
        self.head = nn.Sequential(
            nn.Linear(latent_dim, 128),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(128, n_classes),
        )

    def forward(self, z):
        return self.head(z)


class NASDataset(Dataset):
    """Dataset for NAS fine-tuning."""

    def __init__(self, expr, nas_labels, dataset_idx, fib_stage, indices=None):
        self.indices = indices if indices is not None else np.arange(len(expr))
        self.expr = torch.tensor(expr[self.indices], dtype=torch.float32)
        self.nas_labels = torch.tensor(nas_labels[self.indices], dtype=torch.long)
        self.dataset_idx = torch.tensor(dataset_idx[self.indices], dtype=torch.long)
        self.fib_stage = torch.tensor(fib_stage[self.indices], dtype=torch.long)

    def __len__(self):
        return len(self.indices)

    def __getitem__(self, idx):
        return (self.expr[idx], self.nas_labels[idx],
                self.dataset_idx[idx], self.fib_stage[idx])


def ordinal_contrastive_loss_nas(z, nas_scores, base_margin=0.5):
    """Ordinal contrastive loss using NAS scores (0-8)."""
    mask = nas_scores >= 0
    if mask.sum() < 3:
        return None

    z_labeled = z[mask]
    s_labeled = nas_scores[mask]
    n = z_labeled.shape[0]

    stage_vals = s_labeled.cpu().numpy()
    anchor_idx, pos_idx, neg_idx, margins = [], [], [], []

    for a in range(n):
        sa = stage_vals[a]
        pos_cands = np.where(np.abs(stage_vals - sa) <= 1)[0]
        neg_cands = np.where(np.abs(stage_vals - sa) >= 2)[0]
        pos_cands = pos_cands[pos_cands != a]

        if len(pos_cands) == 0 or len(neg_cands) == 0:
            continue

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


def finetune_vae_for_nas(vae_state_dict, expr, nas_labels, dataset_idx, fib_stage,
                         train_idx, val_idx, n_classes=4, n_epochs=50,
                         lr_encoder=1e-5, lr_head=1e-3):
    """Fine-tune pretrained VAE for NAS prediction."""

    # Load pretrained VAE
    model = ConditionalVAE(
        n_genes=N_GENES, latent_dim=64, stage_embed_dim=8,
        max_stages=10, n_datasets=N_DATASETS, dropout=0.3,
    ).to(DEVICE)
    model.load_state_dict(vae_state_dict)

    # Add NAS head
    nas_head = NASHead(latent_dim=64, n_classes=n_classes).to(DEVICE)

    # Differential learning rates
    optimizer = torch.optim.Adam([
        {"params": model.encoder.parameters(), "lr": lr_encoder},
        {"params": model.fc_mu.parameters(), "lr": lr_encoder},
        {"params": model.fc_logvar.parameters(), "lr": lr_encoder},
        {"params": model.decoder.parameters(), "lr": lr_encoder},
        {"params": model.stage_embed.parameters(), "lr": lr_encoder},
        {"params": model.domain_head.parameters(), "lr": lr_encoder},
        {"params": model.grl.parameters(), "lr": lr_encoder},
        {"params": nas_head.parameters(), "lr": lr_head},
    ])

    train_ds = NASDataset(expr, nas_labels, dataset_idx, fib_stage,
                          indices=train_idx)
    val_ds = NASDataset(expr, nas_labels, dataset_idx, fib_stage,
                        indices=val_idx)
    train_loader = DataLoader(train_ds, batch_size=64, shuffle=True,
                              num_workers=2, pin_memory=True)
    val_loader = DataLoader(val_ds, batch_size=64, shuffle=False,
                            num_workers=2, pin_memory=True)

    best_val_loss = float("inf")
    best_model_state = None
    best_head_state = None
    patience_counter = 0
    patience = 15

    for epoch in range(1, n_epochs + 1):
        model.train()
        nas_head.train()
        train_loss = 0.0
        n_batches = 0

        for x_batch, nas_batch, ds_batch, fib_batch in train_loader:
            x_batch = x_batch.to(DEVICE)
            nas_batch = nas_batch.to(DEVICE)
            ds_batch = ds_batch.to(DEVICE)
            fib_batch = fib_batch.to(DEVICE)

            # Forward through VAE encoder
            x_recon, mu, logvar, z, domain_logits = model(x_batch, fib_batch)

            # Reconstruction + KL
            loss_recon = F.mse_loss(x_recon, x_batch)
            loss_kl = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())

            # NAS ordinal contrastive
            loss_contr = ordinal_contrastive_loss_nas(z, nas_batch)
            if loss_contr is None:
                loss_contr = torch.tensor(0.0, device=DEVICE)

            # NAS classification head (only on labeled samples)
            labeled_mask = nas_batch >= 0
            if labeled_mask.sum() > 0:
                nas_logits = nas_head(mu[labeled_mask])
                loss_cls = F.cross_entropy(nas_logits, nas_batch[labeled_mask])
            else:
                loss_cls = torch.tensor(0.0, device=DEVICE)

            # Domain adversarial
            loss_domain = F.cross_entropy(domain_logits, ds_batch)

            total = (loss_recon + 0.5 * loss_kl + 0.3 * loss_contr
                     + 1.0 * loss_cls + 0.1 * loss_domain)

            optimizer.zero_grad()
            total.backward()
            torch.nn.utils.clip_grad_norm_(
                list(model.parameters()) + list(nas_head.parameters()),
                max_norm=5.0,
            )
            optimizer.step()

            train_loss += total.item()
            n_batches += 1

        train_loss /= max(n_batches, 1)

        # Validation
        model.eval()
        nas_head.eval()
        val_loss = 0.0
        val_batches = 0

        with torch.no_grad():
            for x_batch, nas_batch, ds_batch, fib_batch in val_loader:
                x_batch = x_batch.to(DEVICE)
                nas_batch = nas_batch.to(DEVICE)
                fib_batch = fib_batch.to(DEVICE)

                mu, _ = model.encode(x_batch)
                labeled_mask = nas_batch >= 0
                if labeled_mask.sum() > 0:
                    nas_logits = nas_head(mu[labeled_mask])
                    loss = F.cross_entropy(nas_logits, nas_batch[labeled_mask])
                    val_loss += loss.item()
                    val_batches += 1

        val_loss = val_loss / max(val_batches, 1)

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_model_state = {k: v.cpu().clone()
                                for k, v in model.state_dict().items()}
            best_head_state = {k: v.cpu().clone()
                               for k, v in nas_head.state_dict().items()}
            patience_counter = 0
        else:
            patience_counter += 1

        if patience_counter >= patience:
            break

    return best_model_state, best_head_state


@torch.no_grad()
def embed_and_predict(model_state, head_state, expr, fib_stage, dataset_idx,
                      n_classes, batch_size=256):
    """Generate embeddings and NAS predictions from fine-tuned model."""
    model = ConditionalVAE(
        n_genes=N_GENES, latent_dim=64, stage_embed_dim=8,
        max_stages=10, n_datasets=N_DATASETS, dropout=0.3,
    ).to(DEVICE)
    model.load_state_dict(model_state)
    model.eval()

    nas_head = NASHead(latent_dim=64, n_classes=n_classes).to(DEVICE)
    nas_head.load_state_dict(head_state)
    nas_head.eval()

    all_mu = []
    all_logits = []

    ds = NASDataset(expr, np.zeros(len(expr), dtype=np.int64),
                    dataset_idx, fib_stage)
    loader = DataLoader(ds, batch_size=batch_size, shuffle=False)

    for x_batch, _, ds_batch, fib_batch in loader:
        x_batch = x_batch.to(DEVICE)
        mu, _ = model.encode(x_batch)
        logits = nas_head(mu)
        all_mu.append(mu.cpu().numpy())
        all_logits.append(logits.cpu().numpy())

    return np.vstack(all_mu), np.vstack(all_logits)


def run_finetune():
    """Fine-tune fibrosis VAE for NAS, 5-fold NAS LOCO."""
    print("\n" + "=" * 60)
    print("APPROACH B: Fine-Tune Fibrosis VAE for NAS")
    print("=" * 60)

    # Load pretrained weights
    if not os.path.exists(VAE_PATH):
        print(f"  WARNING: {VAE_PATH} not found, skipping Approach B")
        return []

    vae_state = torch.load(VAE_PATH, map_location="cpu", weights_only=True)
    print(f"  Loaded pretrained VAE: {VAE_PATH}")

    results = []

    for target_name, tinfo in TARGETS.items():
        y = tinfo["y"]
        n_classes = tinfo["n_classes"]
        valid = y >= 0

        print(f"\n  Target: {target_name} ({tinfo['description']})")

        for fold_name in unique_nas_folds:
            test_mask = (loco_nas == fold_name) & valid
            train_all_mask = (loco_nas != fold_name) & (loco_nas != "excluded") & \
                             (loco_nas != "NA") & valid

            if test_mask.sum() < 5 or train_all_mask.sum() < 10:
                continue

            # Further split train into train/val (80/20)
            train_all_idx = np.where(train_all_mask)[0]
            np.random.seed(42)
            perm = np.random.permutation(train_all_idx)
            n_train = int(0.8 * len(perm))
            train_idx = perm[:n_train]
            val_idx = perm[n_train:]

            test_idx = np.where(test_mask)[0]

            print(f"    Fold {fold_name}: train={len(train_idx)}, "
                  f"val={len(val_idx)}, test={len(test_idx)}")

            # Fine-tune
            try:
                best_model_state, best_head_state = finetune_vae_for_nas(
                    copy.deepcopy(vae_state), raw_expr, y.astype(np.int64),
                    dataset_idx, fib_stage.astype(np.int64),
                    train_idx, val_idx, n_classes=n_classes,
                )

                # Evaluate on test
                _, all_logits = embed_and_predict(
                    best_model_state, best_head_state, raw_expr,
                    fib_stage.astype(np.int64), dataset_idx, n_classes,
                )

                test_logits = all_logits[test_idx]
                test_probs = torch.softmax(
                    torch.tensor(test_logits), dim=-1
                ).numpy()
                y_pred = test_logits.argmax(axis=1)
                y_test = y[test_idx].astype(int)

                m = compute_metrics(y_test, y_pred, test_probs, n_classes)
                m.update({
                    "approach": "B_finetune",
                    "target": target_name,
                    "model": "VAE_finetune",
                    "fold": fold_name,
                    "n_train": len(train_idx),
                    "n_test": len(test_idx),
                })
                results.append(m)

                primary = tinfo["metric"]
                print(f"      {primary}={m[primary]:.3f}")

            except Exception as e:
                print(f"      FAILED: {e}")

            # Free GPU memory
            torch.cuda.empty_cache()

    return results


# =========================================================================
# Approach B-scratch: Train NAS VAE from scratch (baseline comparison)
# =========================================================================
def run_scratch():
    """Train VAE from scratch for NAS (no fibrosis pretraining)."""
    print("\n" + "=" * 60)
    print("APPROACH B-scratch: Train NAS VAE from Scratch (Baseline)")
    print("=" * 60)

    results = []

    for target_name, tinfo in TARGETS.items():
        y = tinfo["y"]
        n_classes = tinfo["n_classes"]
        valid = y >= 0

        print(f"\n  Target: {target_name} ({tinfo['description']})")

        for fold_name in unique_nas_folds:
            test_mask = (loco_nas == fold_name) & valid
            train_all_mask = (loco_nas != fold_name) & (loco_nas != "excluded") & \
                             (loco_nas != "NA") & valid

            if test_mask.sum() < 5 or train_all_mask.sum() < 10:
                continue

            train_all_idx = np.where(train_all_mask)[0]
            np.random.seed(42)
            perm = np.random.permutation(train_all_idx)
            n_train = int(0.8 * len(perm))
            train_idx = perm[:n_train]
            val_idx = perm[n_train:]
            test_idx = np.where(test_mask)[0]

            print(f"    Fold {fold_name}: train={len(train_idx)}, test={len(test_idx)}")

            try:
                # Initialize from random weights (no pretraining)
                scratch_model = ConditionalVAE(
                    n_genes=N_GENES, latent_dim=64, stage_embed_dim=8,
                    max_stages=10, n_datasets=N_DATASETS, dropout=0.3,
                )
                scratch_state = {k: v.cpu().clone()
                                 for k, v in scratch_model.state_dict().items()}

                best_model_state, best_head_state = finetune_vae_for_nas(
                    scratch_state, raw_expr, y.astype(np.int64),
                    dataset_idx, fib_stage.astype(np.int64),
                    train_idx, val_idx, n_classes=n_classes,
                    lr_encoder=1e-3, lr_head=1e-3,  # Higher lr for scratch
                )

                _, all_logits = embed_and_predict(
                    best_model_state, best_head_state, raw_expr,
                    fib_stage.astype(np.int64), dataset_idx, n_classes,
                )

                test_logits = all_logits[test_idx]
                test_probs = torch.softmax(
                    torch.tensor(test_logits), dim=-1
                ).numpy()
                y_pred = test_logits.argmax(axis=1)
                y_test = y[test_idx].astype(int)

                m = compute_metrics(y_test, y_pred, test_probs, n_classes)
                m.update({
                    "approach": "B_scratch",
                    "target": target_name,
                    "model": "VAE_scratch",
                    "fold": fold_name,
                    "n_train": len(train_idx),
                    "n_test": len(test_idx),
                })
                results.append(m)

                primary = tinfo["metric"]
                print(f"      {primary}={m[primary]:.3f}")

            except Exception as e:
                print(f"      FAILED: {e}")

            torch.cuda.empty_cache()

    return results


# =========================================================================
# Main
# =========================================================================
def main():
    t0 = time.time()
    np.random.seed(42)
    torch.manual_seed(42)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(42)

    # Run all approaches
    results_a = run_linear_probe()
    results_b = run_finetune()
    results_scratch = run_scratch()

    all_results = results_a + results_b + results_scratch

    # -----------------------------------------------------------------------
    # Save results
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("Saving Results")
    print("=" * 60)

    results_df = pd.DataFrame(all_results)
    results_out = os.path.join(OUTDIR, "fibrosis_transfer_results.csv")
    results_df.to_csv(results_out, index=False)
    print(f"  Results: {results_out} ({len(results_df)} rows)")

    # Build comparison: finetune vs scratch vs linear probe
    comparison_rows = []
    for target_name in TARGETS:
        for approach in ["A_linear_probe", "B_finetune", "B_scratch"]:
            sub = results_df[
                (results_df["target"] == target_name) &
                (results_df["approach"] == approach)
            ]
            if len(sub) == 0:
                continue

            primary = TARGETS[target_name]["metric"]
            comparison_rows.append({
                "target": target_name,
                "approach": approach,
                "primary_metric": primary,
                "mean_value": sub[primary].mean(),
                "std_value": sub[primary].std(),
                "n_folds": len(sub),
                "mean_balanced_accuracy": sub["balanced_accuracy"].mean(),
                "mean_qwk": sub["qwk"].mean(),
                "mean_auroc": sub["auroc"].mean(),
            })

    comparison_df = pd.DataFrame(comparison_rows)
    comparison_out = os.path.join(OUTDIR, "transfer_vs_scratch_comparison.csv")
    comparison_df.to_csv(comparison_out, index=False)
    print(f"  Comparison: {comparison_out} ({len(comparison_df)} rows)")

    # Print summary
    print("\n--- Transfer Learning Summary ---")
    if len(comparison_df) > 0:
        for target_name in TARGETS:
            print(f"\n  {target_name} ({TARGETS[target_name]['metric']}):")
            sub = comparison_df[comparison_df["target"] == target_name]
            for _, row in sub.iterrows():
                print(f"    {row['approach']}: {row['mean_value']:.3f} "
                      f"+/- {row['std_value']:.3f} ({row['n_folds']} folds)")

    dt = time.time() - t0
    print(f"\n{'=' * 60}")
    print(f"93_fibrosis_transfer_nas.py completed in {dt / 60:.1f} min")
    print(f"  Results:    {results_out}")
    print(f"  Comparison: {comparison_out}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
