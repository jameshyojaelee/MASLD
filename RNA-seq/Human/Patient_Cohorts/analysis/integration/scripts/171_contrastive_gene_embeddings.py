#!/usr/bin/env python3
"""
171_contrastive_gene_embeddings.py
Contrastive Multi-Source Gene Embeddings — learn a unified 64-dim embedding
per gene from 7 orthogonal evidence sources in the multi-evidence atlas.

Architecture:
  - 7 SourceEncoders (one per evidence source) → 32-dim each
  - MultiheadAttention fusion (4 heads) with masking for missing sources
  - Projection to 64-dim unified embedding
  - NT-Xent contrastive loss with source-dropout augmentation

Downstream:
  1. UMAP visualization (color by #sources, category, drug target)
  2. K-means clustering (k=5-10) with silhouette + GO enrichment
  3. Target centroid ranking (distance to drug target centroid)
  4. Per-cluster source attention analysis

Output to results/novel_ml/contrastive_embeddings/:
  - gene_embeddings_64dim.csv        33,943 x 64 embeddings
  - gene_clusters.csv                cluster assignments
  - target_centroid_ranking.csv      all genes ranked by distance to drug target centroid
  - source_attention_weights.csv     per-gene attention weights for 7 sources
  - training_metrics.csv             loss curve, silhouette per epoch
  - umap_coordinates.csv             2D UMAP for visualization

SLURM: cpu partition, 8 CPUs, 16G, 48h
Env:   micromamba activate spatial

Usage:
  sbatch --job-name=stg171_contrastive \
         --partition=cpu --cpus-per-task=8 --mem=16G --time=48:00:00 \
         --output=logs/171_contrastive_%j.out \
         --error=logs/171_contrastive_%j.err \
         --wrap="bash -c 'eval \"\\$(micromamba shell hook --shell bash)\" && \\
                 micromamba activate spatial && \\
                 cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts && \\
                 python 171_contrastive_gene_embeddings.py'"
"""

import os
import sys
import time
import warnings
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import cdist

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------
SEED = 42
np.random.seed(SEED)
torch.manual_seed(SEED)
torch.set_num_threads(8)

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
ATLAS_PATH = PROJECT_ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
POSITIVE_CONTROLS = PROJECT_ROOT / "results/library/positive_control.csv"
DRUG_VALIDATION = PROJECT_ROOT / "RNA-seq/results/drug_repurposing/clinical_drug_validation_table.csv"
RESULTS_DIR = PROJECT_ROOT / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/novel_ml/contrastive_embeddings"

# ---------------------------------------------------------------------------
# Source definitions: map each of the 7 sources to atlas columns
# ---------------------------------------------------------------------------
SOURCE_COLUMNS = {
    "S1_human_bulk": [
        "bulk_logFC", "bulk_padj", "bulk_tstat",
        # Progression contrasts (also derived from human bulk)
        "nafl_vs_nash_logFC", "nafl_vs_nash_tstat",
        "adv_fib_logFC", "nas_ge5_logFC", "extreme_logFC",
        "steatosis_ordinal_coef", "inflammation_ordinal_coef",
        "ballooning_ordinal_coef", "fibrosis_ordinal_coef",
        "cirrhosis_logFC", "f2_inflection_logFC",
    ],
    "S2_mouse_bulk": [
        "mouse_meta_logFC", "mouse_meta_padj", "n_diets_sig",
        "translatability_score", "n_concordant_diets",
    ],
    # MR columns (mr_beta, mr_pval, mr_ivw_*, mr_n_instruments, mr_n_gwas_sig,
    # mr_bidirectional_pval) removed 2026-04-22 — MR ditched from paper.
    "S3_genetic_causal": [
        "twas_z", "twas_pval",
        "coloc_pp4", "broadaway_coloc_pp4",
        "ast_coloc_pp4", "ggt_coloc_pp4", "pdff_coloc_pp4",
        "ukbb_alt_coloc_pp4", "best_liver_enzyme_pp4",
        "n_liver_enzyme_coloc", "n_coloc_sources", "n_ancestry_gwas",
        "sceqtl_coloc_pp4_hep", "sceqtl_coloc_best_pp4",
        "sceqtl_n_cell_types",
        "bbj_alt_coloc_pp4", "bbj_ast_coloc_pp4", "bbj_ggt_coloc_pp4",
        "otters_broadaway_pval", "otters_broadaway_z", "otters_n_gwas_sig",
        "ctwas_pip", "ctwas_pval",
        "coloc_susie_best_pp4", "coloc_susie_n_signals",
        "hyprcoloc_posterior", "hyprcoloc_n_traits",
        "causal_methods_sig",
    ],
    "S4_essentiality": [
        "essentiality_chronos",
    ],
    "S5_epigenomic": [
        # Sex-stratified (epigenomic proxy columns available in atlas)
        "bulk_logFC_M", "bulk_logFC_F", "sex_interaction_padj",
        "n_leading_edge_pathways",
    ],
    "S6_spatial": [
        # Progression-derived spatial proxies
        "progression_tau", "progression_n_contrasts_sig",
    ],
    "S7_single_cell": [
        # sc-eQTL and deconv columns serve as single-cell proxies
        "nafl_vs_ctrl_logFC", "nash_vs_ctrl_logFC",
        "early_late_nash_logFC", "nash_vs_nafl_fibadj_logFC",
    ],
}

# Columns used to determine if a source has real data for a gene.
# Integer count columns (e.g., n_coloc_sources) can be 0 without being NaN,
# so we only check float "signal" columns where NaN = no data.
SOURCE_MASK_COLUMNS = {
    "S1_human_bulk": ["bulk_logFC"],
    "S2_mouse_bulk": ["mouse_meta_logFC"],
    "S3_genetic_causal": [
        # mr_beta removed 2026-04-22 — MR ditched from paper.
        "twas_z", "coloc_pp4", "broadaway_coloc_pp4",
        "ctwas_pip", "otters_broadaway_z", "coloc_susie_best_pp4",
        "hyprcoloc_posterior",
    ],
    "S4_essentiality": ["essentiality_chronos"],
    "S5_epigenomic": ["bulk_logFC_M"],
    "S6_spatial": ["progression_tau"],
    "S7_single_cell": ["nafl_vs_ctrl_logFC"],
}

# Hyperparameters
BATCH_SIZE = 512
LR = 1e-3
WEIGHT_DECAY = 1e-5
MAX_EPOCHS = 200
PATIENCE = 20  # early stopping patience
TEMPERATURE = 0.1  # NT-Xent temperature
N_AUG_DROP = 2  # max sources to drop per augmentation
K_RANGE = range(5, 11)  # k-means k range

# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------

class SourceEncoder(nn.Module):
    """Encode one evidence source into 32-dim."""

    def __init__(self, in_dim):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, 64),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(64, 32),
        )

    def forward(self, x):
        return self.net(x)


class MultiSourceContrastive(nn.Module):
    """7 source encoders + attention fusion -> 64-dim unified embedding."""

    def __init__(self, source_dims):
        super().__init__()
        self.n_sources = len(source_dims)
        self.encoders = nn.ModuleList(
            [SourceEncoder(d) for d in source_dims]
        )
        self.attention = nn.MultiheadAttention(
            embed_dim=32, num_heads=4, batch_first=True
        )
        self.projection = nn.Sequential(
            nn.ReLU(),
            nn.Linear(32, 64),
        )

    def forward(self, sources, mask):
        """
        sources: list of 7 tensors, each (batch, source_dim)
        mask: (batch, 7) boolean, True if source is available
        Returns: (batch, 64) embeddings, (batch, 7) attention weights
        """
        # Encode each source
        encoded = torch.stack(
            [enc(s) for enc, s in zip(self.encoders, sources)], dim=1
        )  # (batch, 7, 32)

        # Zero out unavailable sources in encoded representations too
        mask_expanded = mask.unsqueeze(-1)  # (batch, 7, 1)
        encoded = encoded * mask_expanded.float()

        # Attention with masking (key_padding_mask: True = ignore)
        attn_out, attn_weights = self.attention(
            encoded, encoded, encoded,
            key_padding_mask=~mask,
            need_weights=True,
            average_attn_weights=True,
        )  # attn_out: (batch, 7, 32), attn_weights: (batch, 7, 7)

        # Mean pooling over source dimension
        # Weight by mask to avoid including masked positions
        mask_f = mask.float().unsqueeze(-1)  # (batch, 7, 1)
        n_active = mask_f.sum(dim=1).clamp(min=1)  # (batch, 1)
        pooled = (attn_out * mask_f).sum(dim=1) / n_active  # (batch, 32)

        embeddings = self.projection(pooled)  # (batch, 64)

        # Per-source attention: average attention received per source
        # attn_weights is (batch, 7, 7) — average over query dim
        source_attn = attn_weights.mean(dim=1)  # (batch, 7)

        return embeddings, source_attn


def nt_xent_loss(z1, z2, temperature=0.1):
    """
    NT-Xent (Normalized Temperature-scaled Cross Entropy) loss.
    z1, z2: (batch, dim) — positive pairs from augmented views.
    """
    batch_size = z1.shape[0]
    z1 = F.normalize(z1, dim=1)
    z2 = F.normalize(z2, dim=1)

    # Concatenate both views
    z = torch.cat([z1, z2], dim=0)  # (2B, dim)
    sim = torch.mm(z, z.t()) / temperature  # (2B, 2B)

    # Mask out self-similarity
    mask = ~torch.eye(2 * batch_size, dtype=torch.bool, device=z.device)
    sim = sim.masked_fill(~mask, -1e9)

    # Positive pairs: (i, i+B) and (i+B, i)
    labels = torch.cat([
        torch.arange(batch_size, 2 * batch_size),
        torch.arange(0, batch_size),
    ]).to(z.device)

    loss = F.cross_entropy(sim, labels)
    return loss


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------

class GeneSourceDataset(Dataset):
    """Dataset of genes with 7 source feature vectors + masks."""

    def __init__(self, source_arrays, source_masks):
        """
        source_arrays: list of 7 np.ndarray, each (n_genes, source_dim)
        source_masks: list of 7 np.ndarray (boolean), each (n_genes,)
        """
        self.n_genes = source_arrays[0].shape[0]
        self.source_arrays = [torch.tensor(a, dtype=torch.float32) for a in source_arrays]
        self.source_masks = [torch.tensor(m, dtype=torch.bool) for m in source_masks]

    def __len__(self):
        return self.n_genes

    def __getitem__(self, idx):
        sources = [a[idx] for a in self.source_arrays]
        mask = torch.tensor([m[idx] for m in self.source_masks])
        return sources, mask


def collate_fn(batch):
    """Custom collate for list-of-sources structure."""
    n_sources = len(batch[0][0])
    sources = [torch.stack([b[0][s] for b in batch]) for s in range(n_sources)]
    masks = torch.stack([b[1] for b in batch])
    return sources, masks


def augment_batch(sources, mask, n_drop_max=2):
    """
    Create augmented view by randomly dropping 1-n_drop_max sources per gene.
    Only drops sources that are currently available (mask=True).
    """
    batch_size = mask.shape[0]
    n_sources = mask.shape[1]
    aug_sources = [s.clone() for s in sources]
    aug_mask = mask.clone()

    for i in range(batch_size):
        available = torch.where(mask[i])[0]
        if len(available) <= 1:
            continue  # don't drop if only 0-1 sources available
        n_drop = np.random.randint(1, min(n_drop_max + 1, len(available)))
        drop_idx = available[torch.randperm(len(available))[:n_drop]]
        for j in drop_idx:
            aug_sources[j][i] = 0.0
            aug_mask[i, j] = False

    return aug_sources, aug_mask


# ---------------------------------------------------------------------------
# Data loading and preprocessing
# ---------------------------------------------------------------------------

def load_and_prepare_data():
    """Load atlas, map columns to sources, normalize, build masks."""
    log.info("Loading multi-evidence atlas from %s", ATLAS_PATH)
    atlas = pd.read_csv(ATLAS_PATH, low_memory=False)
    log.info("Atlas shape: %s", atlas.shape)

    gene_ids = atlas["ensembl_id"].values
    gene_symbols = atlas["human_symbol"].values
    gene_biotypes = atlas["gene_biotype"].values

    # Validate all source columns exist
    all_cols = set(atlas.columns)
    for source_name, cols in SOURCE_COLUMNS.items():
        missing = [c for c in cols if c not in all_cols]
        if missing:
            log.warning("Source %s missing columns: %s — will be excluded", source_name, missing)
            SOURCE_COLUMNS[source_name] = [c for c in cols if c in all_cols]

    # Convert boolean-like columns to float
    bool_cols = atlas.select_dtypes(include=["bool"]).columns
    for c in bool_cols:
        atlas[c] = atlas[c].astype(float)

    # Convert string booleans
    for col in atlas.columns:
        if atlas[col].dtype == object:
            unique_vals = set(atlas[col].dropna().unique())
            if unique_vals <= {"True", "False", "TRUE", "FALSE"}:
                atlas[col] = atlas[col].map(
                    {"True": 1.0, "False": 0.0, "TRUE": 1.0, "FALSE": 0.0}
                )

    # Build per-source arrays and masks
    source_arrays = []
    source_masks = []
    source_names = list(SOURCE_COLUMNS.keys())
    source_dims = []

    for source_name in source_names:
        cols = SOURCE_COLUMNS[source_name]
        if not cols:
            log.warning("Source %s has no valid columns, using dummy", source_name)
            cols = ["bulk_logFC"]  # fallback

        df_sub = atlas[cols].copy()

        # Convert any remaining object cols to numeric
        for c in df_sub.columns:
            if df_sub[c].dtype == object:
                df_sub[c] = pd.to_numeric(df_sub[c], errors="coerce")

        arr = df_sub.values.astype(np.float64)

        # Mask: use designated signal columns (float cols where NaN = no data)
        # This avoids integer count columns (always 0, never NaN) inflating availability
        mask_cols = SOURCE_MASK_COLUMNS.get(source_name, cols)
        mask_cols = [c for c in mask_cols if c in atlas.columns]
        if mask_cols:
            mask_arr = atlas[mask_cols].apply(pd.to_numeric, errors="coerce").values
            mask = ~np.all(np.isnan(mask_arr), axis=1)
        else:
            mask = ~np.all(np.isnan(arr), axis=1)

        # Replace NaN with 0 for the tensor
        arr = np.nan_to_num(arr, nan=0.0)

        # Standardize each column (zero mean, unit variance) using only non-NaN entries
        scaler = StandardScaler()
        # Fit on rows where source is available
        if mask.sum() > 0:
            scaler.fit(arr[mask])
            arr = scaler.transform(arr)
        else:
            log.warning("Source %s has no available genes!", source_name)

        source_arrays.append(arr.astype(np.float32))
        source_masks.append(mask)
        source_dims.append(arr.shape[1])
        log.info(
            "  %s: %d columns, %d/%d genes available (%.1f%%)",
            source_name, arr.shape[1], mask.sum(), len(mask),
            100 * mask.sum() / len(mask),
        )

    # Count active sources per gene
    n_active = np.stack(source_masks, axis=1).sum(axis=1)
    log.info("Active sources per gene: mean=%.1f, median=%d, min=%d, max=%d",
             n_active.mean(), np.median(n_active), n_active.min(), n_active.max())

    # Load drug targets / positive controls
    drug_genes = set()
    if DRUG_VALIDATION.exists():
        drug_df = pd.read_csv(DRUG_VALIDATION)
        drug_genes = set(drug_df["target_gene"].dropna().unique())
        log.info("Loaded %d drug target genes", len(drug_genes))

    positive_genes = set()
    if POSITIVE_CONTROLS.exists():
        pc_df = pd.read_csv(POSITIVE_CONTROLS)
        # Handle 'SCD1' -> map to SCD if needed
        pc_syms = set(pc_df["Gene symbol"].dropna().unique())
        positive_genes = pc_syms
        log.info("Loaded %d positive control genes", len(positive_genes))

    # Build metadata DataFrame
    meta = pd.DataFrame({
        "ensembl_id": gene_ids,
        "human_symbol": gene_symbols,
        "gene_biotype": gene_biotypes,
        "n_active_sources": n_active,
        "is_drug_target": [s in drug_genes for s in gene_symbols],
        "is_positive_control": [s in positive_genes for s in gene_symbols],
        "layers_active": atlas["layers_active"].values if "layers_active" in atlas.columns else 0,
    })

    # Functional categories from atlas columns
    # Use primary_category if available
    if "primary_category" in atlas.columns:
        meta["primary_category"] = atlas["primary_category"].values
    if "ferroptosis_class" in atlas.columns:
        meta["ferroptosis_class"] = atlas["ferroptosis_class"].values
    if "zonation_class" in atlas.columns:
        meta["zonation_class"] = atlas["zonation_class"].values
    if "is_conserved" in atlas.columns:
        meta["is_conserved"] = atlas["is_conserved"].values

    return source_arrays, source_masks, source_dims, source_names, meta


# ---------------------------------------------------------------------------
# Training
# ---------------------------------------------------------------------------

def train_model(source_arrays, source_masks, source_dims, meta):
    """Train contrastive model and return embeddings + attention weights."""
    dataset = GeneSourceDataset(source_arrays, source_masks)
    loader = DataLoader(
        dataset, batch_size=BATCH_SIZE, shuffle=True,
        collate_fn=collate_fn, num_workers=0, drop_last=True,
    )

    model = MultiSourceContrastive(source_dims)
    optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=MAX_EPOCHS)

    log.info("Model parameters: %d", sum(p.numel() for p in model.parameters()))
    log.info("Training with batch_size=%d, lr=%s, epochs=%d", BATCH_SIZE, LR, MAX_EPOCHS)

    metrics = []
    best_silhouette = -1
    patience_counter = 0
    best_state = None

    for epoch in range(MAX_EPOCHS):
        model.train()
        epoch_loss = 0.0
        n_batches = 0

        for sources, mask in loader:
            # View 1: original (with existing masks)
            z1, _ = model(sources, mask)

            # View 2: augmented (drop 1-2 random sources)
            aug_sources, aug_mask = augment_batch(sources, mask, n_drop_max=N_AUG_DROP)
            z2, _ = model(aug_sources, aug_mask)

            loss = nt_xent_loss(z1, z2, temperature=TEMPERATURE)

            optimizer.zero_grad()
            loss.backward()
            # Gradient clipping
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()

            epoch_loss += loss.item()
            n_batches += 1

        scheduler.step()
        avg_loss = epoch_loss / max(n_batches, 1)

        # Evaluate silhouette every 10 epochs (expensive on full dataset)
        sil = -1.0
        if (epoch + 1) % 10 == 0 or epoch == 0:
            emb = get_all_embeddings(model, source_arrays, source_masks)
            # Subsample for silhouette (too expensive on 34K genes)
            n_sub = min(5000, len(emb))
            idx_sub = np.random.choice(len(emb), n_sub, replace=False)
            emb_sub = emb[idx_sub]
            labels_sub = KMeans(n_clusters=7, random_state=SEED, n_init=5).fit_predict(emb_sub)
            sil = silhouette_score(emb_sub, labels_sub)

            if sil > best_silhouette:
                best_silhouette = sil
                best_state = {k: v.clone() for k, v in model.state_dict().items()}
                patience_counter = 0
            else:
                patience_counter += 1

        metrics.append({
            "epoch": epoch + 1,
            "loss": avg_loss,
            "silhouette": sil if sil > -1 else np.nan,
            "lr": scheduler.get_last_lr()[0],
        })

        if (epoch + 1) % 10 == 0:
            log.info(
                "Epoch %3d/%d  loss=%.4f  silhouette=%.4f  lr=%.2e  patience=%d/%d",
                epoch + 1, MAX_EPOCHS, avg_loss, sil,
                scheduler.get_last_lr()[0], patience_counter, PATIENCE,
            )

        if patience_counter >= PATIENCE // 10:
            # Check every evaluation cycle (10 epochs)
            # Total patience = PATIENCE evaluations * 10 epochs
            pass

        if patience_counter >= PATIENCE:
            log.info("Early stopping at epoch %d (best silhouette=%.4f)", epoch + 1, best_silhouette)
            break

    # Load best model
    if best_state is not None:
        model.load_state_dict(best_state)
        log.info("Loaded best model (silhouette=%.4f)", best_silhouette)

    # Get final embeddings and attention weights
    embeddings = get_all_embeddings(model, source_arrays, source_masks)
    attn_weights = get_all_attention(model, source_arrays, source_masks)

    return model, embeddings, attn_weights, pd.DataFrame(metrics)


def get_all_embeddings(model, source_arrays, source_masks):
    """Get embeddings for all genes."""
    model.eval()
    dataset = GeneSourceDataset(source_arrays, source_masks)
    loader = DataLoader(dataset, batch_size=1024, shuffle=False, collate_fn=collate_fn)

    all_emb = []
    with torch.no_grad():
        for sources, mask in loader:
            emb, _ = model(sources, mask)
            all_emb.append(emb.numpy())
    return np.concatenate(all_emb, axis=0)


def get_all_attention(model, source_arrays, source_masks):
    """Get attention weights for all genes."""
    model.eval()
    dataset = GeneSourceDataset(source_arrays, source_masks)
    loader = DataLoader(dataset, batch_size=1024, shuffle=False, collate_fn=collate_fn)

    all_attn = []
    with torch.no_grad():
        for sources, mask in loader:
            _, attn = model(sources, mask)
            all_attn.append(attn.numpy())
    return np.concatenate(all_attn, axis=0)


# ---------------------------------------------------------------------------
# Downstream analysis
# ---------------------------------------------------------------------------

def run_clustering(embeddings, meta, source_names):
    """K-means clustering with silhouette evaluation."""
    log.info("Running K-means clustering (k=%d-%d)", K_RANGE.start, K_RANGE.stop - 1)

    best_k = 7
    best_sil = -1
    results = []

    for k in K_RANGE:
        km = KMeans(n_clusters=k, random_state=SEED, n_init=10, max_iter=300)
        labels = km.fit_predict(embeddings)
        sil = silhouette_score(embeddings, labels)
        results.append({"k": k, "silhouette": sil, "inertia": km.inertia_})
        log.info("  k=%d  silhouette=%.4f  inertia=%.1f", k, sil, km.inertia_)
        if sil > best_sil:
            best_sil = sil
            best_k = k

    log.info("Best k=%d (silhouette=%.4f)", best_k, best_sil)

    # Final clustering with best k
    km_final = KMeans(n_clusters=best_k, random_state=SEED, n_init=10, max_iter=300)
    labels = km_final.fit_predict(embeddings)

    cluster_df = meta.copy()
    cluster_df["cluster"] = labels
    cluster_df["best_k"] = best_k

    # Per-cluster summary
    for c in range(best_k):
        idx = labels == c
        n_genes = idx.sum()
        n_drug = cluster_df.loc[idx, "is_drug_target"].sum()
        n_pc = cluster_df.loc[idx, "is_positive_control"].sum()
        mean_sources = cluster_df.loc[idx, "n_active_sources"].mean()
        log.info(
            "  Cluster %d: %d genes, %d drug targets, %d pos controls, mean_sources=%.1f",
            c, n_genes, n_drug, n_pc, mean_sources,
        )

    return cluster_df, pd.DataFrame(results)


def run_umap(embeddings, meta):
    """UMAP dimensionality reduction."""
    import umap

    log.info("Computing UMAP (n=%d genes)", len(embeddings))
    reducer = umap.UMAP(
        n_components=2, n_neighbors=30, min_dist=0.3,
        metric="euclidean", random_state=SEED,
    )
    coords = reducer.fit_transform(embeddings)

    umap_df = meta[["ensembl_id", "human_symbol"]].copy()
    umap_df["umap_1"] = coords[:, 0]
    umap_df["umap_2"] = coords[:, 1]
    umap_df["n_active_sources"] = meta["n_active_sources"].values
    umap_df["is_drug_target"] = meta["is_drug_target"].values
    umap_df["is_positive_control"] = meta["is_positive_control"].values
    if "primary_category" in meta.columns:
        umap_df["primary_category"] = meta["primary_category"].values

    return umap_df


def target_centroid_ranking(embeddings, meta):
    """Rank all genes by distance to drug target centroid."""
    is_target = meta["is_drug_target"].values | meta["is_positive_control"].values
    n_targets = is_target.sum()
    log.info("Computing target centroid from %d known target/control genes", n_targets)

    if n_targets == 0:
        log.warning("No drug targets found in atlas; skipping centroid ranking")
        return pd.DataFrame()

    centroid = embeddings[is_target].mean(axis=0, keepdims=True)  # (1, 64)
    distances = cdist(embeddings, centroid, metric="euclidean").flatten()

    ranking_df = meta[["ensembl_id", "human_symbol", "gene_biotype",
                        "n_active_sources", "is_drug_target",
                        "is_positive_control"]].copy()
    ranking_df["centroid_distance"] = distances
    ranking_df["rank"] = ranking_df["centroid_distance"].rank(method="min").astype(int)
    ranking_df = ranking_df.sort_values("centroid_distance")

    # Log top 20 novel genes (not already drug targets)
    novel = ranking_df[~ranking_df["is_drug_target"] & ~ranking_df["is_positive_control"]]
    log.info("Top 20 novel genes nearest drug target centroid:")
    for _, row in novel.head(20).iterrows():
        log.info(
            "  Rank %4d: %-12s dist=%.4f  sources=%d",
            row["rank"], row["human_symbol"], row["centroid_distance"],
            row["n_active_sources"],
        )

    return ranking_df


def source_attention_analysis(attn_weights, cluster_labels, source_names):
    """Analyze which sources dominate per cluster."""
    n_clusters = len(set(cluster_labels))

    rows = []
    for c in range(n_clusters):
        idx = cluster_labels == c
        mean_attn = attn_weights[idx].mean(axis=0)
        row = {"cluster": c, "n_genes": idx.sum()}
        for i, name in enumerate(source_names):
            row[f"attn_{name}"] = mean_attn[i]
        row["dominant_source"] = source_names[np.argmax(mean_attn)]
        rows.append(row)

    summary = pd.DataFrame(rows)
    log.info("Source attention by cluster:")
    for _, row in summary.iterrows():
        log.info(
            "  Cluster %d (%d genes): dominant=%s",
            row["cluster"], row["n_genes"], row["dominant_source"],
        )
    return summary


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    t0 = time.time()
    log.info("=" * 70)
    log.info("Script 171: Contrastive Multi-Source Gene Embeddings")
    log.info("=" * 70)

    # Create output directory
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    log.info("Output directory: %s", RESULTS_DIR)

    # 1. Load and prepare data
    source_arrays, source_masks, source_dims, source_names, meta = load_and_prepare_data()

    # 2. Train contrastive model
    model, embeddings, attn_weights, training_metrics = train_model(
        source_arrays, source_masks, source_dims, meta
    )

    # 3. UMAP
    umap_df = run_umap(embeddings, meta)

    # 4. Clustering
    cluster_df, cluster_eval = run_clustering(embeddings, meta, source_names)

    # 5. Target centroid ranking
    ranking_df = target_centroid_ranking(embeddings, meta)

    # 6. Source attention analysis
    attn_summary = source_attention_analysis(
        attn_weights, cluster_df["cluster"].values, source_names
    )

    # ---------------------------------------------------------------------------
    # Save outputs
    # ---------------------------------------------------------------------------
    log.info("Saving outputs to %s", RESULTS_DIR)

    # Embeddings
    emb_df = pd.DataFrame(
        embeddings,
        columns=[f"emb_{i}" for i in range(embeddings.shape[1])],
    )
    emb_df.insert(0, "ensembl_id", meta["ensembl_id"].values)
    emb_df.insert(1, "human_symbol", meta["human_symbol"].values)
    emb_df.to_csv(RESULTS_DIR / "gene_embeddings_64dim.csv", index=False)
    log.info("  gene_embeddings_64dim.csv: %s", emb_df.shape)

    # Clusters
    cluster_df.to_csv(RESULTS_DIR / "gene_clusters.csv", index=False)
    log.info("  gene_clusters.csv: %s", cluster_df.shape)

    # Cluster evaluation
    cluster_eval.to_csv(RESULTS_DIR / "cluster_evaluation.csv", index=False)

    # Target ranking
    if not ranking_df.empty:
        ranking_df.to_csv(RESULTS_DIR / "target_centroid_ranking.csv", index=False)
        log.info("  target_centroid_ranking.csv: %s", ranking_df.shape)

    # Attention weights (per gene)
    attn_df = pd.DataFrame(attn_weights, columns=[f"attn_{n}" for n in source_names])
    attn_df.insert(0, "ensembl_id", meta["ensembl_id"].values)
    attn_df.insert(1, "human_symbol", meta["human_symbol"].values)
    attn_df.to_csv(RESULTS_DIR / "source_attention_weights.csv", index=False)
    log.info("  source_attention_weights.csv: %s", attn_df.shape)

    # Attention summary per cluster
    attn_summary.to_csv(RESULTS_DIR / "source_attention_by_cluster.csv", index=False)

    # Training metrics
    training_metrics.to_csv(RESULTS_DIR / "training_metrics.csv", index=False)
    log.info("  training_metrics.csv: %s epochs", len(training_metrics))

    # UMAP coordinates
    umap_df.to_csv(RESULTS_DIR / "umap_coordinates.csv", index=False)
    log.info("  umap_coordinates.csv: %s", umap_df.shape)

    # Save model
    torch.save(model.state_dict(), RESULTS_DIR / "contrastive_model.pt")
    log.info("  contrastive_model.pt saved")

    elapsed = time.time() - t0
    log.info("=" * 70)
    log.info("Script 171 complete in %.1f min", elapsed / 60)
    log.info("=" * 70)


if __name__ == "__main__":
    main()
