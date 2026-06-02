#!/usr/bin/env python3
"""
71_attention_panel.py
Transformer-based gene panel selection via CLS-token attention.

Architecture:
  Gene tokenizer: Each of 3000 genes -> Linear(1->128) + learned gene positional embedding
  Transformer encoder: 2 layers, 4 heads, dim=128, ff_dim=256, dropout=0.1
  CLS token: prepended, used for classification
  Classification heads: Linear(128->n_classes) per target

Attention extraction:
  Average CLS->gene attention across heads and layers -> per-gene score
  Rank genes -> attention panels (10, 20, 50 genes)
  5 random seeds for stability assessment

Input:  results/staging_classifier/prepared_data.h5
Output: results/staging_classifier/
  - attention_panel_10.csv, attention_panel_20.csv, attention_panel_50.csv
  - attention_weights_all.csv (per-gene, mean/std across seeds)
  - attention_stability.csv (CV per gene)
  - attention_vs_stability_comparison.csv
  - transformer_loco_results.csv

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
from collections import defaultdict

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
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

print("=" * 60)
print("71: Transformer Attention-Based Gene Panel Selection")
print("=" * 60)
print(f"Device: {DEVICE}")
print(f"Input:  {H5_PATH}")
print()

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
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

N_SAMPLES, N_GENES = expr.shape
print(f"  Samples: {N_SAMPLES}, Genes: {N_GENES}")

# Targets
fib_stage = meta["fib_stage"].astype(np.int64)
nas_group4 = meta["nas_group4"].astype(np.int64)
loco_folds_fib = meta["loco_fold_fibrosis"]
loco_folds_nas = meta["loco_fold_nas"]

unique_folds_fib = sorted(set(f for f in loco_folds_fib if f not in ("excluded", "NA")))
unique_folds_nas = sorted(set(f for f in loco_folds_nas if f not in ("excluded", "NA")))

print(f"  Fibrosis-labeled: {(fib_stage >= 0).sum()}")
print(f"  NAS-labeled: {(nas_group4 >= 0).sum()}")
print(f"  LOCO fib folds: {len(unique_folds_fib)}")
print(f"  LOCO NAS folds: {len(unique_folds_nas)}")

# Model hyperparams
SEEDS = [42, 123, 456, 789, 1024]
N_HEADS = 4
N_LAYERS = 2
D_MODEL = 128
D_FF = 256
DROPOUT = 0.1
BATCH_SIZE = 32
N_EPOCHS = 100
PATIENCE = 15
LR = 1e-4


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------
class GeneExprDataset(Dataset):
    """Dataset for gene-tokenized transformer."""

    def __init__(self, expr, fib_stage, nas_group4, indices=None):
        self.indices = indices if indices is not None else np.arange(len(expr))
        self.expr = torch.tensor(expr[self.indices], dtype=torch.float32)
        self.fib = torch.tensor(fib_stage[self.indices], dtype=torch.long)
        self.nas = torch.tensor(nas_group4[self.indices], dtype=torch.long)

    def __len__(self):
        return len(self.indices)

    def __getitem__(self, idx):
        return self.expr[idx], self.fib[idx], self.nas[idx]


# ---------------------------------------------------------------------------
# Transformer Model
# ---------------------------------------------------------------------------
class GeneTransformer(nn.Module):
    """
    Transformer with gene tokenization and CLS token for classification.

    Each gene is a "token": expression value * linear projection + gene identity embedding.
    CLS token is prepended. Sequence length = N_GENES + 1.
    """

    def __init__(self, n_genes, d_model=128, n_heads=4, n_layers=2,
                 d_ff=256, dropout=0.1, n_fib_classes=5, n_nas_classes=4):
        super().__init__()
        self.n_genes = n_genes
        self.d_model = d_model

        # Gene value projection: scalar -> d_model
        self.value_proj = nn.Linear(1, d_model)

        # Gene identity embeddings (not positional -- based on gene identity)
        self.gene_embed = nn.Embedding(n_genes, d_model)

        # CLS token (learnable)
        self.cls_token = nn.Parameter(torch.randn(1, 1, d_model) * 0.02)

        # Transformer encoder
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model, nhead=n_heads, dim_feedforward=d_ff,
            dropout=dropout, batch_first=True, activation="relu",
        )
        self.transformer = nn.TransformerEncoder(
            encoder_layer, num_layers=n_layers,
        )

        # Classification heads
        self.head_fib = nn.Linear(d_model, n_fib_classes)
        self.head_nas = nn.Linear(d_model, n_nas_classes)

        # Layer norm before classification
        self.ln = nn.LayerNorm(d_model)

    def forward(self, x, return_attention=False):
        """
        x: (batch, n_genes) expression values
        Returns: fib_logits, nas_logits, [attention_weights if requested]
        """
        B = x.shape[0]

        # Project each gene's expression value
        x_proj = self.value_proj(x.unsqueeze(-1))  # (B, n_genes, d_model)

        # Add gene identity embeddings
        gene_ids = torch.arange(self.n_genes, device=x.device)
        x_proj = x_proj + self.gene_embed(gene_ids).unsqueeze(0)  # broadcast

        # Prepend CLS token
        cls = self.cls_token.expand(B, -1, -1)  # (B, 1, d_model)
        seq = torch.cat([cls, x_proj], dim=1)  # (B, n_genes+1, d_model)

        # Transformer encoding
        if return_attention:
            # Extract attention weights manually by hooking into layers
            attn_weights = []

            def hook_fn(module, input, output):
                # For TransformerEncoderLayer, we need to get attention from
                # the self_attn sub-module. We register on self_attn directly.
                pass

            # Use a custom forward to capture attention
            encoded = seq
            for layer in self.transformer.layers:
                # Manual forward to extract attention
                # self_attn expects (query, key, value)
                src = encoded
                src2, attn_w = layer.self_attn(
                    src, src, src, need_weights=True, average_attn_weights=False,
                )
                attn_weights.append(attn_w.detach())  # (B, n_heads, seq_len, seq_len)
                src = src + layer.dropout1(src2)
                src = layer.norm1(src)
                src2 = layer.linear2(layer.dropout(layer.activation(layer.linear1(src))))
                src = src + layer.dropout2(src2)
                src = layer.norm2(src)
                encoded = src
        else:
            encoded = self.transformer(seq)
            attn_weights = None

        # CLS token output
        cls_out = self.ln(encoded[:, 0, :])  # (B, d_model)

        fib_logits = self.head_fib(cls_out)
        nas_logits = self.head_nas(cls_out)

        if return_attention:
            return fib_logits, nas_logits, attn_weights
        return fib_logits, nas_logits


# ---------------------------------------------------------------------------
# Class weights
# ---------------------------------------------------------------------------
def compute_class_weights(y, n_classes):
    """Compute balanced class weights."""
    valid = y >= 0
    if valid.sum() == 0:
        return torch.ones(n_classes)
    y_valid = y[valid]
    classes, counts = np.unique(y_valid, return_counts=True)
    n = len(y_valid)
    weights = torch.ones(n_classes)
    for cls, cnt in zip(classes, counts):
        if int(cls) < n_classes:
            weights[int(cls)] = n / (n_classes * cnt)
    return weights


# ---------------------------------------------------------------------------
# Training one model instance
# ---------------------------------------------------------------------------
def train_transformer(model, train_loader, val_loader, fib_weights, nas_weights,
                      n_epochs=100, patience=15, lr=1e-4):
    """Train transformer with early stopping on validation loss."""
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    fib_weights = fib_weights.to(DEVICE)
    nas_weights = nas_weights.to(DEVICE)

    best_val_loss = float("inf")
    best_state = None
    patience_counter = 0
    history = []

    for epoch in range(1, n_epochs + 1):
        model.train()
        train_loss = 0.0
        n_batches = 0

        for x_batch, fib_batch, nas_batch in train_loader:
            x_batch = x_batch.to(DEVICE)
            fib_batch = fib_batch.to(DEVICE)
            nas_batch = nas_batch.to(DEVICE)

            fib_logits, nas_logits = model(x_batch)

            loss = torch.tensor(0.0, device=DEVICE)
            n_targets = 0

            # Fibrosis loss (only labeled samples)
            fib_mask = fib_batch >= 0
            if fib_mask.sum() > 0:
                loss = loss + F.cross_entropy(
                    fib_logits[fib_mask], fib_batch[fib_mask], weight=fib_weights,
                )
                n_targets += 1

            # NAS loss (only labeled samples)
            nas_mask = nas_batch >= 0
            if nas_mask.sum() > 0:
                loss = loss + F.cross_entropy(
                    nas_logits[nas_mask], nas_batch[nas_mask], weight=nas_weights,
                )
                n_targets += 1

            if n_targets > 0:
                loss = loss / n_targets
                optimizer.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                optimizer.step()

            train_loss += loss.item()
            n_batches += 1

        train_loss /= max(n_batches, 1)

        # Validation
        model.eval()
        val_loss = 0.0
        val_batches = 0
        with torch.no_grad():
            for x_batch, fib_batch, nas_batch in val_loader:
                x_batch = x_batch.to(DEVICE)
                fib_batch = fib_batch.to(DEVICE)
                nas_batch = nas_batch.to(DEVICE)

                fib_logits, nas_logits = model(x_batch)
                loss = torch.tensor(0.0, device=DEVICE)
                n_t = 0
                fib_mask = fib_batch >= 0
                if fib_mask.sum() > 0:
                    loss = loss + F.cross_entropy(
                        fib_logits[fib_mask], fib_batch[fib_mask], weight=fib_weights,
                    )
                    n_t += 1
                nas_mask = nas_batch >= 0
                if nas_mask.sum() > 0:
                    loss = loss + F.cross_entropy(
                        nas_logits[nas_mask], nas_batch[nas_mask], weight=nas_weights,
                    )
                    n_t += 1
                if n_t > 0:
                    loss = loss / n_t
                val_loss += loss.item()
                val_batches += 1

        val_loss /= max(val_batches, 1)
        history.append({"epoch": epoch, "train_loss": train_loss, "val_loss": val_loss})

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            patience_counter = 0
        else:
            patience_counter += 1

        if patience_counter >= patience:
            break

    return best_state, history


# ---------------------------------------------------------------------------
# Extract attention weights from CLS token
# ---------------------------------------------------------------------------
@torch.no_grad()
def extract_attention(model, expr_tensor, batch_size=32):
    """
    Extract CLS->gene attention weights averaged across heads and layers.
    Returns: (n_samples, n_genes) attention matrix.
    """
    model.eval()
    ds = torch.utils.data.TensorDataset(expr_tensor)
    loader = DataLoader(ds, batch_size=batch_size, shuffle=False)

    all_attn = []
    for (x_batch,) in loader:
        x_batch = x_batch.to(DEVICE)
        B = x_batch.shape[0]

        # Forward with attention capture
        _, _, attn_weights = model(x_batch, return_attention=True)

        # attn_weights: list of (B, n_heads, seq_len, seq_len) per layer
        # We want CLS (position 0) -> gene tokens (positions 1:)
        batch_attn = []
        for layer_attn in attn_weights:
            # (B, n_heads, seq_len, seq_len) -> CLS row -> gene columns
            cls_attn = layer_attn[:, :, 0, 1:]  # (B, n_heads, n_genes)
            batch_attn.append(cls_attn)

        # Average across heads and layers
        stacked = torch.stack(batch_attn, dim=0)  # (n_layers, B, n_heads, n_genes)
        avg_attn = stacked.mean(dim=(0, 2))  # (B, n_genes)
        all_attn.append(avg_attn.cpu().numpy())

    return np.vstack(all_attn)


# ---------------------------------------------------------------------------
# Evaluate LOCO fold
# ---------------------------------------------------------------------------
def evaluate_loco_fold(model, X_test, fib_test, nas_test, batch_size=64):
    """Evaluate trained model on held-out fold."""
    from sklearn.metrics import balanced_accuracy_score, cohen_kappa_score

    model.eval()
    ds = torch.utils.data.TensorDataset(
        torch.tensor(X_test, dtype=torch.float32),
        torch.tensor(fib_test, dtype=torch.long),
        torch.tensor(nas_test, dtype=torch.long),
    )
    loader = DataLoader(ds, batch_size=batch_size, shuffle=False)

    fib_preds, fib_trues = [], []
    nas_preds, nas_trues = [], []

    with torch.no_grad():
        for x_b, f_b, n_b in loader:
            x_b = x_b.to(DEVICE)
            fib_logits, nas_logits = model(x_b)

            fib_mask = f_b >= 0
            if fib_mask.sum() > 0:
                fib_preds.extend(fib_logits[fib_mask].argmax(dim=1).cpu().numpy())
                fib_trues.extend(f_b[fib_mask].numpy())

            nas_mask = n_b >= 0
            if nas_mask.sum() > 0:
                nas_preds.extend(nas_logits[nas_mask].argmax(dim=1).cpu().numpy())
                nas_trues.extend(n_b[nas_mask].numpy())

    metrics = {}
    if len(fib_trues) >= 2:
        metrics["fib_balanced_acc"] = balanced_accuracy_score(fib_trues, fib_preds)
        metrics["fib_kappa"] = cohen_kappa_score(fib_trues, fib_preds, weights="quadratic")
        metrics["fib_n"] = len(fib_trues)
    if len(nas_trues) >= 2:
        metrics["nas_balanced_acc"] = balanced_accuracy_score(nas_trues, nas_preds)
        metrics["nas_kappa"] = cohen_kappa_score(nas_trues, nas_preds, weights="quadratic")
        metrics["nas_n"] = len(nas_trues)

    return metrics


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()

    # Storage for per-seed attention weights
    seed_attention_maps = {}  # seed -> (n_genes,) averaged attention scores
    all_loco_results = []

    for seed in SEEDS:
        print(f"\n{'='*60}")
        print(f"  Seed: {seed}")
        print(f"{'='*60}")

        torch.manual_seed(seed)
        np.random.seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)

        # Class weights (computed from all labeled samples)
        fib_weights = compute_class_weights(fib_stage, 5)
        nas_weights = compute_class_weights(nas_group4, 4)

        # ----- LOCO cross-validation -----
        for fold_name in unique_folds_fib:
            fold_test_mask = np.array([f == fold_name for f in loco_folds_fib])
            fold_train_all = ~fold_test_mask

            # Further split train -> train/val (80/20)
            train_cands = np.where(fold_train_all)[0]
            rng = np.random.RandomState(seed)
            perm = rng.permutation(train_cands)
            n_tr = int(0.8 * len(perm))
            train_idx = perm[:n_tr]
            val_idx = perm[n_tr:]

            train_ds = GeneExprDataset(expr, fib_stage, nas_group4, indices=train_idx)
            val_ds = GeneExprDataset(expr, fib_stage, nas_group4, indices=val_idx)
            train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True,
                                      num_workers=2, pin_memory=True, drop_last=False)
            val_loader = DataLoader(val_ds, batch_size=BATCH_SIZE, shuffle=False,
                                    num_workers=2, pin_memory=True)

            model = GeneTransformer(
                n_genes=N_GENES, d_model=D_MODEL, n_heads=N_HEADS,
                n_layers=N_LAYERS, d_ff=D_FF, dropout=DROPOUT,
            ).to(DEVICE)

            best_state, hist = train_transformer(
                model, train_loader, val_loader, fib_weights, nas_weights,
                n_epochs=N_EPOCHS, patience=PATIENCE, lr=LR,
            )

            model.load_state_dict(best_state)
            model.to(DEVICE)

            # Evaluate on held-out fold
            test_idx = np.where(fold_test_mask)[0]
            metrics = evaluate_loco_fold(
                model, expr[test_idx], fib_stage[test_idx], nas_group4[test_idx],
            )
            metrics["seed"] = seed
            metrics["fold"] = fold_name
            metrics["n_epochs"] = len(hist)
            all_loco_results.append(metrics)

            print(f"    Fold {fold_name}: epochs={len(hist)}", end="")
            if "fib_balanced_acc" in metrics:
                print(f"  fib_bacc={metrics['fib_balanced_acc']:.3f} "
                      f"kappa={metrics['fib_kappa']:.3f}", end="")
            if "nas_balanced_acc" in metrics:
                print(f"  nas_bacc={metrics['nas_balanced_acc']:.3f}", end="")
            print()

        # ----- Extract attention from a model trained on ALL data -----
        print(f"\n  Training full model for attention extraction (seed {seed})...")

        # 80/20 split for full model
        all_idx = np.arange(N_SAMPLES)
        rng = np.random.RandomState(seed)
        perm = rng.permutation(all_idx)
        n_tr = int(0.8 * len(perm))

        full_train_ds = GeneExprDataset(expr, fib_stage, nas_group4, indices=perm[:n_tr])
        full_val_ds = GeneExprDataset(expr, fib_stage, nas_group4, indices=perm[n_tr:])
        full_train_loader = DataLoader(full_train_ds, batch_size=BATCH_SIZE, shuffle=True,
                                        num_workers=2, pin_memory=True)
        full_val_loader = DataLoader(full_val_ds, batch_size=BATCH_SIZE, shuffle=False,
                                      num_workers=2, pin_memory=True)

        full_model = GeneTransformer(
            n_genes=N_GENES, d_model=D_MODEL, n_heads=N_HEADS,
            n_layers=N_LAYERS, d_ff=D_FF, dropout=DROPOUT,
        ).to(DEVICE)

        best_full_state, _ = train_transformer(
            full_model, full_train_loader, full_val_loader,
            fib_weights, nas_weights,
            n_epochs=N_EPOCHS, patience=PATIENCE, lr=LR,
        )
        full_model.load_state_dict(best_full_state)
        full_model.to(DEVICE)

        # Extract attention from labeled samples only (fibrosis-staged)
        labeled_mask = fib_stage >= 0
        labeled_expr = torch.tensor(expr[labeled_mask], dtype=torch.float32)
        attn_map = extract_attention(full_model, labeled_expr, batch_size=32)
        # Average across samples -> per-gene score
        gene_attn = attn_map.mean(axis=0)  # (n_genes,)
        seed_attention_maps[seed] = gene_attn

        print(f"  Attention extracted: {gene_attn.shape} (max={gene_attn.max():.4f})")

        # Free GPU memory
        del model, full_model
        torch.cuda.empty_cache()

    # -----------------------------------------------------------------------
    # Aggregate attention across seeds
    # -----------------------------------------------------------------------
    print("\n--- Aggregating attention across seeds ---")

    attn_matrix = np.stack([seed_attention_maps[s] for s in SEEDS], axis=0)  # (n_seeds, n_genes)
    attn_mean = attn_matrix.mean(axis=0)
    attn_std = attn_matrix.std(axis=0)
    attn_cv = np.where(attn_mean > 0, attn_std / attn_mean, 0.0)

    # Build attention summary DataFrame
    attn_df = pd.DataFrame({
        "gene": gene_names,
        "attention_mean": attn_mean,
        "attention_std": attn_std,
        "attention_cv": attn_cv,
    })
    for i, seed in enumerate(SEEDS):
        attn_df[f"seed_{seed}"] = attn_matrix[i]

    attn_df = attn_df.sort_values("attention_mean", ascending=False).reset_index(drop=True)
    attn_df["rank"] = range(1, len(attn_df) + 1)

    attn_path = os.path.join(OUTDIR, "attention_weights_all.csv")
    attn_df.to_csv(attn_path, index=False)
    print(f"  Attention weights saved: {attn_path}")

    # -----------------------------------------------------------------------
    # Gene panels
    # -----------------------------------------------------------------------
    for panel_size in [10, 20, 50]:
        panel = attn_df.head(panel_size)[["gene", "attention_mean", "attention_std",
                                            "attention_cv", "rank"]].copy()
        panel_path = os.path.join(OUTDIR, f"attention_panel_{panel_size}.csv")
        panel.to_csv(panel_path, index=False)
        genes_str = ", ".join(panel["gene"].tolist()[:10])
        if panel_size > 10:
            genes_str += ", ..."
        print(f"  Panel {panel_size}: {panel_path}")
        print(f"    Top genes: {genes_str}")

    # -----------------------------------------------------------------------
    # Stability analysis
    # -----------------------------------------------------------------------
    print("\n--- Stability Analysis ---")
    stability_df = attn_df[["gene", "attention_mean", "attention_cv", "rank"]].copy()
    stability_df["is_stable"] = stability_df["attention_cv"] < 0.3  # CV < 30%
    stability_df["is_unstable"] = stability_df["attention_cv"] > 0.5  # CV > 50%

    stability_path = os.path.join(OUTDIR, "attention_stability.csv")
    stability_df.to_csv(stability_path, index=False)
    print(f"  Stability saved: {stability_path}")
    print(f"    Stable genes (CV<0.3): {stability_df['is_stable'].sum()}")
    print(f"    Unstable genes (CV>0.5): {stability_df['is_unstable'].sum()}")

    # -----------------------------------------------------------------------
    # Compare attention-selected vs stability-selected panels
    # -----------------------------------------------------------------------
    print("\n--- Attention vs Stability Panel Comparison ---")

    # Attention panel: top N by mean attention
    # Stability panel: top N among stable genes (CV<0.3), ranked by attention
    comparison_rows = []
    for panel_size in [10, 20, 50]:
        attn_panel = set(attn_df.head(panel_size)["gene"])

        stable_ranked = stability_df[stability_df["is_stable"]].sort_values(
            "attention_mean", ascending=False
        )
        stab_panel = set(stable_ranked.head(panel_size)["gene"])

        overlap = len(attn_panel & stab_panel)
        jaccard = overlap / len(attn_panel | stab_panel) if len(attn_panel | stab_panel) > 0 else 0

        comparison_rows.append({
            "panel_size": panel_size,
            "attention_only": len(attn_panel - stab_panel),
            "stability_only": len(stab_panel - attn_panel),
            "overlap": overlap,
            "jaccard": jaccard,
            "attention_genes": ";".join(sorted(attn_panel)),
            "stability_genes": ";".join(sorted(stab_panel)),
        })
        print(f"    Panel {panel_size}: overlap={overlap}, Jaccard={jaccard:.3f}")

    comp_df = pd.DataFrame(comparison_rows)
    comp_path = os.path.join(OUTDIR, "attention_vs_stability_comparison.csv")
    comp_df.to_csv(comp_path, index=False)

    # -----------------------------------------------------------------------
    # LOCO results
    # -----------------------------------------------------------------------
    loco_df = pd.DataFrame(all_loco_results)
    loco_path = os.path.join(OUTDIR, "transformer_loco_results.csv")
    loco_df.to_csv(loco_path, index=False)
    print(f"\nLOCO results saved: {loco_path} ({len(loco_df)} rows)")

    # Summary statistics
    if "fib_balanced_acc" in loco_df.columns:
        print("\n  Fibrosis balanced accuracy across seeds:")
        for seed in SEEDS:
            seed_rows = loco_df[loco_df["seed"] == seed]
            fib_vals = seed_rows["fib_balanced_acc"].dropna()
            if len(fib_vals) > 0:
                print(f"    Seed {seed}: mean={fib_vals.mean():.3f} "
                      f"std={fib_vals.std():.3f} ({len(fib_vals)} folds)")
        all_fib = loco_df["fib_balanced_acc"].dropna()
        print(f"    Overall: mean={all_fib.mean():.3f} std={all_fib.std():.3f}")

    if "nas_balanced_acc" in loco_df.columns:
        all_nas = loco_df["nas_balanced_acc"].dropna()
        if len(all_nas) > 0:
            print(f"\n  NAS balanced accuracy: mean={all_nas.mean():.3f} "
                  f"std={all_nas.std():.3f}")

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    dt = time.time() - t0
    print(f"\n{'=' * 60}")
    print(f"71_attention_panel.py completed in {dt / 60:.1f} min")
    print(f"  Panels:     {OUTDIR}/attention_panel_{{10,20,50}}.csv")
    print(f"  Weights:    {attn_path}")
    print(f"  Stability:  {stability_path}")
    print(f"  Comparison: {comp_path}")
    print(f"  LOCO:       {loco_path}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
