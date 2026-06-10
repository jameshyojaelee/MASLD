#!/usr/bin/env python3
"""
171v2_contrastive_real_sources.py
Contrastive Multi-Source Gene Embeddings v2 — ONLY 4 truly independent sources.

Fixes over v1 (Script 171):
  - REMOVED S5/S6/S7 (epigenomic, spatial, single-cell) — all were derivatives
    of S1 human bulk, not independent evidence layers.
  - 4 source encoders (S1 Human bulk, S2 Mouse bulk, S3 Genetic causal,
    S4 Essentiality) — no fake sources.
  - Embedding dim: 32 (not 64) — appropriate for 4 sparse sources.
  - NaN handling: mask=False when ALL source columns are NaN; NaN->0 only in
    the tensor after masking is set. No imputation leakage.
  - Temperature: 0.5 (not 0.1) — less aggressive contrastive loss for sparse data.
  - Early stopping: properly tracks silhouette every 5 epochs and stops.
  - Batch size: 256 (not 512).
  - Baselines: PCA + K-means, layers_active correlation.
  - GO enrichment: Fisher's exact test with MSigDB Hallmark gene sets per cluster.

Architecture:
  - 4 SourceEncoders (one per evidence source) -> 32-dim each
  - MultiheadAttention fusion (4 heads) with masking for missing sources
  - Projection to 32-dim unified embedding
  - NT-Xent contrastive loss with source-dropout augmentation

Output to results/novel_ml/contrastive_embeddings_v2/:
  - gene_embeddings_32dim.csv        33,943 x 32 embeddings
  - gene_clusters.csv                cluster assignments + metadata
  - target_centroid_ranking.csv      all genes ranked by distance to drug target centroid
  - source_attention_weights.csv     per-gene attention weights for 4 sources
  - training_metrics.csv             loss curve, silhouette per epoch
  - umap_coordinates.csv             2D UMAP for visualization
  - cluster_evaluation.csv           K-means silhouette per k
  - source_attention_by_cluster.csv  per-cluster attention summary
  - baseline_comparison.csv          PCA baseline vs contrastive silhouette
  - layers_active_correlation.csv    layers_active vs centroid distance
  - cluster_go_enrichment.csv        Fisher's test Hallmark enrichment per cluster
  - contrastive_model.pt             trained model weights

SLURM: cpu partition, 8 CPUs, 16G, 48h
Env:   micromamba activate spatial

Usage:
  sbatch --job-name=stg171v2_contrastive \
         --partition=cpu --cpus-per-task=8 --mem=16G --time=48:00:00 \
         --output=logs/171v2_contrastive_%j.out \
         --error=logs/171v2_contrastive_%j.err \
         --wrap="bash -c 'eval \"\\$(micromamba shell hook --shell bash)\" && \\
                 micromamba activate spatial && \\
                 cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts && \\
                 python 171v2_contrastive_real_sources.py'"
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
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import cdist
from scipy.stats import fisher_exact, spearmanr

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
DRUG_VALIDATION = (
    PROJECT_ROOT / "RNA-seq/results/drug_repurposing/clinical_drug_validation_table.csv"
)
GMT_PATH = (
    PROJECT_ROOT
    / "Analysis/downstream_analysis/pathway_analysis/data/msigdb.v2025.1.Hs.symbols.gmt"
)
RESULTS_DIR = (
    PROJECT_ROOT
    / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/novel_ml/contrastive_embeddings_v2"
)

# ---------------------------------------------------------------------------
# Source definitions: ONLY 4 truly independent evidence sources
# ---------------------------------------------------------------------------
# S3 columns: all numeric columns that are genotype/GWAS/eQTL-derived.
# String/categorical S3 columns are one-hot encoded below.
S3_NUMERIC_COLUMNS = [
    # MR (Mendelian randomization) columns removed 2026-04-22 — MR ditched from paper.
    # TWAS
    "twas_z",
    "twas_pval",
    # COLOC (colocalization across GWAS sources)
    "coloc_pp4",
    "broadaway_coloc_pp4",
    "ast_coloc_pp4",
    "ggt_coloc_pp4",
    "pdff_coloc_pp4",
    "ukbb_alt_coloc_pp4",
    "finngen_nafld_coloc_pp4",
    "finngen_nash_coloc_pp4",
    "finngen_hcc_coloc_pp4",
    "bbj_alt_coloc_pp4",
    "bbj_ast_coloc_pp4",
    "bbj_ggt_coloc_pp4",
    "ghouse_hcc_coloc_pp4",
    "panukbb_afr_alt_coloc_pp4",
    "panukbb_afr_ast_coloc_pp4",
    "panukbb_afr_ggt_coloc_pp4",
    "panukbb_csa_alt_coloc_pp4",
    "panukbb_csa_ast_coloc_pp4",
    "panukbb_csa_ggt_coloc_pp4",
    # sc-eQTL COLOC
    "sceqtl_coloc_pp4_hep",
    "sceqtl_coloc_best_pp4",
    "sceqtl_n_cell_types",
    # SuSiE COLOC
    "coloc_susie_best_pp4",
    "coloc_susie_n_signals",
    # HyPrColoc
    "hyprcoloc_posterior",
    "hyprcoloc_n_traits",
    # Summary counts
    "n_liver_enzyme_coloc",
    "n_coloc_sources",
    "n_ancestry_gwas",
    # ieQTL
    "ieqtl_interaction_pval",
    # Pleiotropy
    "pleiotropy_n_traits",
    "pleiotropy_n_domains",
    # Enhanced MR columns removed 2026-04-22 — MR ditched from paper.
    # OTTERS TWAS
    "otters_broadaway_pval",
    "otters_broadaway_z",
    "otters_n_gwas_sig",
    # cTWAS
    "ctwas_pip",
    "ctwas_pval",
    # Causal summary
    "causal_methods_sig",
    # Progression TWAS (GWAS-based, independent of S1 dream)
    "progression_twas_n_contrasts",
]

# Boolean S3 columns (will be cast to 0/1 float)
S3_BOOLEAN_COLUMNS = [
    "cross_ancestry_coloc_replication",
    "zenodo_nafld_coloc",
    "ieqtl_disease_interaction",
    "is_nafld_specific",
]

SOURCE_COLUMNS = {
    "S1_human_bulk": [
        "bulk_logFC",
        "bulk_padj",
        "bulk_tstat",
    ],
    "S2_mouse_bulk": [
        "mouse_meta_logFC",
        "mouse_meta_padj",
        "n_diets_sig",
    ],
    "S3_genetic_causal": None,  # built dynamically below
    "S4_essentiality": [
        "essentiality_chronos",
    ],
}

# Columns used for masking: a source is "available" for a gene if
# at least one of these signal columns is non-NaN.
SOURCE_MASK_COLUMNS = {
    "S1_human_bulk": ["bulk_logFC"],
    "S2_mouse_bulk": ["mouse_meta_logFC"],
    "S3_genetic_causal": [
        "mr_beta",
        "twas_z",
        "coloc_pp4",
        "broadaway_coloc_pp4",
        "ctwas_pip",
        "otters_broadaway_z",
        "coloc_susie_best_pp4",
        "hyprcoloc_posterior",
        "ieqtl_interaction_pval",
        "sceqtl_coloc_pp4_hep",
    ],
    "S4_essentiality": ["essentiality_chronos"],
}

# Hyperparameters (v2 fixes)
EMBED_DIM = 32
BATCH_SIZE = 256
LR = 1e-3
WEIGHT_DECAY = 1e-5
MAX_EPOCHS = 200
PATIENCE = 20  # early stopping: 20 evaluation cycles * 5 epochs = 100 epochs max wait
EVAL_EVERY = 5  # evaluate silhouette every N epochs
TEMPERATURE = 0.5
N_AUG_DROP = 2  # max sources to drop per augmentation
K_RANGE = range(4, 11)  # k-means k range (4 sources -> start at k=4)


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------


class SourceEncoder(nn.Module):
    """Encode one evidence source into 32-dim."""

    def __init__(self, in_dim: int):
        super().__init__()
        hidden = max(32, in_dim)  # at least 32 hidden units
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(hidden, 32),
        )

    def forward(self, x):
        return self.net(x)


class MultiSourceContrastive(nn.Module):
    """4 source encoders + attention fusion -> 32-dim unified embedding."""

    def __init__(self, source_dims: list[int]):
        super().__init__()
        self.n_sources = len(source_dims)
        self.encoders = nn.ModuleList(
            [SourceEncoder(d) for d in source_dims]
        )
        self.attention = nn.MultiheadAttention(
            embed_dim=32, num_heads=4, batch_first=True
        )
        self.projection = nn.Sequential(
            nn.LayerNorm(32),
            nn.ReLU(),
            nn.Linear(32, EMBED_DIM),
        )

    def forward(self, sources: list[torch.Tensor], mask: torch.Tensor):
        """
        sources: list of N_SOURCES tensors, each (batch, source_dim)
        mask: (batch, N_SOURCES) boolean, True if source is available
        Returns: (batch, EMBED_DIM) embeddings, (batch, N_SOURCES) attention weights
        """
        # Encode each source
        encoded = torch.stack(
            [enc(s) for enc, s in zip(self.encoders, sources)], dim=1
        )  # (batch, N_SOURCES, 32)

        # Zero out unavailable sources
        mask_expanded = mask.unsqueeze(-1).float()  # (batch, N_SOURCES, 1)
        encoded = encoded * mask_expanded

        # Attention with masking (key_padding_mask: True = ignore)
        attn_out, attn_weights = self.attention(
            encoded,
            encoded,
            encoded,
            key_padding_mask=~mask,
            need_weights=True,
            average_attn_weights=True,
        )  # attn_out: (batch, N_SOURCES, 32), attn_weights: (batch, N_SOURCES, N_SOURCES)

        # Masked mean pooling over source dimension
        mask_f = mask.float().unsqueeze(-1)  # (batch, N_SOURCES, 1)
        n_active = mask_f.sum(dim=1).clamp(min=1)  # (batch, 1)
        pooled = (attn_out * mask_f).sum(dim=1) / n_active  # (batch, 32)

        embeddings = self.projection(pooled)  # (batch, EMBED_DIM)

        # Per-source attention: average attention received per source
        source_attn = attn_weights.mean(dim=1)  # (batch, N_SOURCES)

        return embeddings, source_attn


def nt_xent_loss(z1: torch.Tensor, z2: torch.Tensor, temperature: float = 0.5):
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
    eye_mask = torch.eye(2 * batch_size, dtype=torch.bool, device=z.device)
    sim = sim.masked_fill(eye_mask, -1e9)

    # Positive pairs: (i, i+B) and (i+B, i)
    labels = torch.cat(
        [
            torch.arange(batch_size, 2 * batch_size),
            torch.arange(0, batch_size),
        ]
    ).to(z.device)

    loss = F.cross_entropy(sim, labels)
    return loss


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------


class GeneSourceDataset(Dataset):
    """Dataset of genes with 4 source feature vectors + masks."""

    def __init__(self, source_arrays: list[np.ndarray], source_masks: list[np.ndarray]):
        self.n_genes = source_arrays[0].shape[0]
        self.source_arrays = [
            torch.tensor(a, dtype=torch.float32) for a in source_arrays
        ]
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


def augment_batch(
    sources: list[torch.Tensor], mask: torch.Tensor, n_drop_max: int = 2
):
    """
    Create augmented view by randomly dropping 1-n_drop_max sources per gene.
    Only drops sources that are currently available (mask=True).
    Guarantees at least 1 source remains.
    """
    batch_size = mask.shape[0]
    aug_sources = [s.clone() for s in sources]
    aug_mask = mask.clone()

    for i in range(batch_size):
        available = torch.where(mask[i])[0]
        if len(available) <= 1:
            continue  # keep at least 1 source
        max_drop = min(n_drop_max, len(available) - 1)
        n_drop = np.random.randint(1, max_drop + 1)
        drop_idx = available[torch.randperm(len(available))[:n_drop]]
        for j in drop_idx:
            aug_sources[j.item()][i] = 0.0
            aug_mask[i, j.item()] = False

    return aug_sources, aug_mask


# ---------------------------------------------------------------------------
# GMT parsing for GO enrichment
# ---------------------------------------------------------------------------


def load_hallmark_gmt(gmt_path: Path) -> dict[str, set[str]]:
    """Load Hallmark gene sets from MSigDB GMT file."""
    gene_sets = {}
    with open(gmt_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            name = parts[0]
            if name.startswith("HALLMARK_"):
                genes = set(parts[2:])  # skip name and URL
                gene_sets[name] = genes
    log.info("Loaded %d Hallmark gene sets from GMT", len(gene_sets))
    return gene_sets


def fisher_enrichment(
    cluster_genes: set[str],
    background_genes: set[str],
    gene_sets: dict[str, set[str]],
    min_overlap: int = 3,
) -> pd.DataFrame:
    """Fisher's exact test enrichment for a cluster against Hallmark gene sets."""
    N = len(background_genes)
    n_cluster = len(cluster_genes)
    rows = []

    for gs_name, gs_genes in gene_sets.items():
        gs_in_bg = gs_genes & background_genes
        K = len(gs_in_bg)
        if K < 5:
            continue

        overlap = cluster_genes & gs_in_bg
        k = len(overlap)
        if k < min_overlap:
            continue

        # 2x2 contingency table
        a = k  # in cluster AND in gene set
        b = n_cluster - k  # in cluster, NOT in gene set
        c = K - k  # NOT in cluster, in gene set
        d = N - n_cluster - K + k  # neither

        _, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
        odds_ratio = (a * d) / max(b * c, 1)

        rows.append(
            {
                "gene_set": gs_name,
                "overlap": k,
                "gene_set_size": K,
                "cluster_size": n_cluster,
                "odds_ratio": odds_ratio,
                "pval": pval,
                "overlap_genes": ";".join(sorted(overlap)),
            }
        )

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows).sort_values("pval")

    # BH FDR correction
    n_tests = len(df)
    df["rank"] = range(1, n_tests + 1)
    df["padj"] = df["pval"] * n_tests / df["rank"]
    df["padj"] = df["padj"].clip(upper=1.0)
    # Ensure monotonicity
    df["padj"] = df["padj"].iloc[::-1].cummin().iloc[::-1]
    df = df.drop(columns=["rank"])

    return df


# ---------------------------------------------------------------------------
# Data loading and preprocessing
# ---------------------------------------------------------------------------


def load_and_prepare_data():
    """Load atlas, map columns to 4 independent sources, normalize, build masks."""
    log.info("Loading multi-evidence atlas from %s", ATLAS_PATH)
    atlas = pd.read_csv(ATLAS_PATH, low_memory=False)
    log.info("Atlas shape: %s", atlas.shape)

    gene_ids = atlas["ensembl_id"].values
    gene_symbols = atlas["human_symbol"].values
    gene_biotypes = atlas["gene_biotype"].values if "gene_biotype" in atlas.columns else np.array(["unknown"] * len(atlas))

    # -------------------------------------------------------------------
    # Build S3 columns dynamically: numeric + boolean columns
    # -------------------------------------------------------------------
    s3_cols = []

    # Numeric columns
    for c in S3_NUMERIC_COLUMNS:
        if c in atlas.columns:
            s3_cols.append(c)
        else:
            log.warning("S3 numeric column not found: %s", c)

    # Boolean columns -> cast to float
    for c in S3_BOOLEAN_COLUMNS:
        if c in atlas.columns:
            # Cast bool/object to float
            col_name = f"{c}_float"
            if atlas[c].dtype == bool:
                atlas[col_name] = atlas[c].astype(float)
            elif atlas[c].dtype == object:
                atlas[col_name] = atlas[c].map(
                    {"True": 1.0, "False": 0.0, "TRUE": 1.0, "FALSE": 0.0, True: 1.0, False: 0.0}
                )
            else:
                atlas[col_name] = pd.to_numeric(atlas[c], errors="coerce")
            s3_cols.append(col_name)

    log.info("S3 genetic causal: %d columns", len(s3_cols))
    SOURCE_COLUMNS["S3_genetic_causal"] = s3_cols

    # -------------------------------------------------------------------
    # Build per-source arrays and masks
    # -------------------------------------------------------------------
    source_arrays = []
    source_masks = []
    source_names = list(SOURCE_COLUMNS.keys())
    source_dims = []

    for source_name in source_names:
        cols = SOURCE_COLUMNS[source_name]
        if not cols:
            raise ValueError(f"Source {source_name} has no columns!")

        # Extract and convert to numeric
        df_sub = atlas[cols].copy()
        for c in df_sub.columns:
            if df_sub[c].dtype == object:
                df_sub[c] = pd.to_numeric(df_sub[c], errors="coerce")
            elif df_sub[c].dtype == bool:
                df_sub[c] = df_sub[c].astype(float)

        arr = df_sub.values.astype(np.float64)

        # --- MASK: source is available if at least one designated signal
        #     column is non-NaN. This is computed BEFORE NaN->0 imputation. ---
        mask_col_names = SOURCE_MASK_COLUMNS.get(source_name, cols)
        # Filter to columns actually present
        mask_col_names = [c for c in mask_col_names if c in atlas.columns]
        if mask_col_names:
            mask_arr = atlas[mask_col_names].apply(pd.to_numeric, errors="coerce").values
            gene_mask = ~np.all(np.isnan(mask_arr), axis=1)
        else:
            gene_mask = ~np.all(np.isnan(arr), axis=1)

        # --- Standardize ONLY on genes where source is available ---
        # Replace NaN with 0 AFTER computing mask
        arr_clean = np.nan_to_num(arr, nan=0.0)

        scaler = StandardScaler()
        if gene_mask.sum() > 0:
            scaler.fit(arr_clean[gene_mask])
            arr_scaled = scaler.transform(arr_clean)
        else:
            log.warning("Source %s has no available genes!", source_name)
            arr_scaled = arr_clean

        # For genes where source is NOT available, zero out the features.
        # This ensures the encoder gets zeros + mask=False, so the attention
        # layer ignores them.
        arr_scaled[~gene_mask] = 0.0

        source_arrays.append(arr_scaled.astype(np.float32))
        source_masks.append(gene_mask)
        source_dims.append(arr_scaled.shape[1])
        log.info(
            "  %s: %d columns, %d/%d genes available (%.1f%%)",
            source_name,
            arr_scaled.shape[1],
            gene_mask.sum(),
            len(gene_mask),
            100 * gene_mask.sum() / len(gene_mask),
        )

    # Count active sources per gene
    n_active = np.stack(source_masks, axis=1).astype(int).sum(axis=1)
    log.info(
        "Active sources per gene: mean=%.2f, median=%d, min=%d, max=%d",
        n_active.mean(),
        int(np.median(n_active)),
        n_active.min(),
        n_active.max(),
    )

    # Distribution of active source counts
    for k in range(5):
        count = (n_active == k).sum()
        log.info("  %d sources: %d genes (%.1f%%)", k, count, 100 * count / len(n_active))

    # -------------------------------------------------------------------
    # Load drug targets / positive controls
    # -------------------------------------------------------------------
    drug_genes = set()
    if DRUG_VALIDATION.exists():
        drug_df = pd.read_csv(DRUG_VALIDATION)
        if "target_gene" in drug_df.columns:
            drug_genes = set(drug_df["target_gene"].dropna().unique())
        log.info("Loaded %d drug target genes", len(drug_genes))

    positive_genes = set()
    if POSITIVE_CONTROLS.exists():
        pc_df = pd.read_csv(POSITIVE_CONTROLS)
        if "Gene symbol" in pc_df.columns:
            positive_genes = set(pc_df["Gene symbol"].dropna().unique())
        log.info("Loaded %d positive control genes", len(positive_genes))

    # -------------------------------------------------------------------
    # Build metadata DataFrame
    # -------------------------------------------------------------------
    meta = pd.DataFrame(
        {
            "ensembl_id": gene_ids,
            "human_symbol": gene_symbols,
            "gene_biotype": gene_biotypes,
            "n_active_sources": n_active,
            "is_drug_target": [s in drug_genes for s in gene_symbols],
            "is_positive_control": [s in positive_genes for s in gene_symbols],
            "layers_active": (
                atlas["layers_active"].values
                if "layers_active" in atlas.columns
                else np.zeros(len(atlas))
            ),
        }
    )

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


def get_all_embeddings(
    model: nn.Module,
    source_arrays: list[np.ndarray],
    source_masks: list[np.ndarray],
) -> np.ndarray:
    """Get embeddings for all genes (inference mode)."""
    model.eval()
    dataset = GeneSourceDataset(source_arrays, source_masks)
    loader = DataLoader(dataset, batch_size=1024, shuffle=False, collate_fn=collate_fn)

    all_emb = []
    with torch.no_grad():
        for sources, mask in loader:
            emb, _ = model(sources, mask)
            all_emb.append(emb.numpy())
    return np.concatenate(all_emb, axis=0)


def get_all_attention(
    model: nn.Module,
    source_arrays: list[np.ndarray],
    source_masks: list[np.ndarray],
) -> np.ndarray:
    """Get per-gene attention weights (inference mode)."""
    model.eval()
    dataset = GeneSourceDataset(source_arrays, source_masks)
    loader = DataLoader(dataset, batch_size=1024, shuffle=False, collate_fn=collate_fn)

    all_attn = []
    with torch.no_grad():
        for sources, mask in loader:
            _, attn = model(sources, mask)
            all_attn.append(attn.numpy())
    return np.concatenate(all_attn, axis=0)


def evaluate_silhouette(
    embeddings: np.ndarray, k: int = 6, n_sub: int = 5000
) -> float:
    """Compute silhouette score on a subsample."""
    n = len(embeddings)
    n_sub = min(n_sub, n)
    idx = np.random.choice(n, n_sub, replace=False)
    emb_sub = embeddings[idx]
    labels = KMeans(n_clusters=k, random_state=SEED, n_init=5).fit_predict(emb_sub)
    return silhouette_score(emb_sub, labels)


def train_model(
    source_arrays: list[np.ndarray],
    source_masks: list[np.ndarray],
    source_dims: list[int],
    meta: pd.DataFrame,
):
    """Train contrastive model and return embeddings + attention weights."""
    dataset = GeneSourceDataset(source_arrays, source_masks)
    loader = DataLoader(
        dataset,
        batch_size=BATCH_SIZE,
        shuffle=True,
        collate_fn=collate_fn,
        num_workers=0,
        drop_last=True,
    )

    model = MultiSourceContrastive(source_dims)
    optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=MAX_EPOCHS)

    n_params = sum(p.numel() for p in model.parameters())
    log.info("Model parameters: %d", n_params)
    log.info(
        "Training: batch_size=%d, lr=%s, temperature=%.2f, embed_dim=%d, epochs=%d",
        BATCH_SIZE, LR, TEMPERATURE, EMBED_DIM, MAX_EPOCHS,
    )

    metrics = []
    best_silhouette = -1.0
    patience_counter = 0
    best_state = None

    for epoch in range(MAX_EPOCHS):
        model.train()
        epoch_loss = 0.0
        n_batches = 0

        for sources, mask in loader:
            # View 1: original
            z1, _ = model(sources, mask)

            # View 2: augmented (drop 1-2 random sources)
            aug_sources, aug_mask = augment_batch(sources, mask, n_drop_max=N_AUG_DROP)
            z2, _ = model(aug_sources, aug_mask)

            loss = nt_xent_loss(z1, z2, temperature=TEMPERATURE)

            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()

            epoch_loss += loss.item()
            n_batches += 1

        scheduler.step()
        avg_loss = epoch_loss / max(n_batches, 1)

        # Evaluate silhouette every EVAL_EVERY epochs
        sil = np.nan
        if (epoch + 1) % EVAL_EVERY == 0 or epoch == 0:
            emb = get_all_embeddings(model, source_arrays, source_masks)
            sil = evaluate_silhouette(emb, k=6)

            if sil > best_silhouette:
                best_silhouette = sil
                best_state = {k: v.clone() for k, v in model.state_dict().items()}
                patience_counter = 0
            else:
                patience_counter += 1

            log.info(
                "Epoch %3d/%d  loss=%.4f  silhouette=%.4f  lr=%.2e  patience=%d/%d",
                epoch + 1,
                MAX_EPOCHS,
                avg_loss,
                sil,
                scheduler.get_last_lr()[0],
                patience_counter,
                PATIENCE,
            )

            # --- REAL early stopping ---
            if patience_counter >= PATIENCE:
                log.info(
                    "Early stopping at epoch %d (best silhouette=%.4f)",
                    epoch + 1,
                    best_silhouette,
                )
                break
        else:
            if (epoch + 1) % 20 == 0:
                log.info(
                    "Epoch %3d/%d  loss=%.4f  lr=%.2e",
                    epoch + 1, MAX_EPOCHS, avg_loss, scheduler.get_last_lr()[0],
                )

        metrics.append(
            {
                "epoch": epoch + 1,
                "loss": avg_loss,
                "silhouette": sil,
                "lr": scheduler.get_last_lr()[0],
            }
        )

    # Load best model
    if best_state is not None:
        model.load_state_dict(best_state)
        log.info("Loaded best model (silhouette=%.4f)", best_silhouette)
    else:
        log.warning("No best state saved; using final model")

    # Get final embeddings and attention weights
    embeddings = get_all_embeddings(model, source_arrays, source_masks)
    attn_weights = get_all_attention(model, source_arrays, source_masks)

    return model, embeddings, attn_weights, pd.DataFrame(metrics)


# ---------------------------------------------------------------------------
# Baselines
# ---------------------------------------------------------------------------


def run_pca_baseline(
    source_arrays: list[np.ndarray], source_masks: list[np.ndarray]
) -> tuple[np.ndarray, float, pd.DataFrame]:
    """
    PCA baseline: concatenate all raw (scaled) features -> PCA to EMBED_DIM -> K-means.
    Returns PCA embeddings, best silhouette, and evaluation DataFrame.
    """
    log.info("=" * 40)
    log.info("BASELINE: PCA on concatenated raw features")

    # Concatenate all source arrays
    X = np.concatenate(source_arrays, axis=1)
    log.info("  Concatenated shape: %s", X.shape)

    pca = PCA(n_components=EMBED_DIM, random_state=SEED)
    X_pca = pca.fit_transform(X)
    log.info(
        "  PCA explained variance: %.3f (cumulative)",
        pca.explained_variance_ratio_.sum(),
    )

    # K-means evaluation
    best_k_pca = 6
    best_sil_pca = -1.0
    pca_eval_rows = []

    for k in K_RANGE:
        km = KMeans(n_clusters=k, random_state=SEED, n_init=10)
        labels = km.fit_predict(X_pca)
        sil = silhouette_score(X_pca, labels)
        pca_eval_rows.append({"method": "PCA", "k": k, "silhouette": sil})
        if sil > best_sil_pca:
            best_sil_pca = sil
            best_k_pca = k

    log.info("  PCA best k=%d, silhouette=%.4f", best_k_pca, best_sil_pca)
    return X_pca, best_sil_pca, pd.DataFrame(pca_eval_rows)


def run_layers_active_baseline(
    embeddings: np.ndarray, meta: pd.DataFrame
) -> pd.DataFrame:
    """
    Baseline: correlate simple `layers_active` count with embedding distance
    to drug target centroid.
    """
    log.info("=" * 40)
    log.info("BASELINE: layers_active vs centroid distance")

    is_target = meta["is_drug_target"].values | meta["is_positive_control"].values
    n_targets = is_target.sum()

    if n_targets == 0:
        log.warning("No drug targets; skipping layers_active baseline")
        return pd.DataFrame()

    centroid = embeddings[is_target].mean(axis=0, keepdims=True)
    distances = cdist(embeddings, centroid, metric="euclidean").flatten()

    layers = meta["layers_active"].values.astype(float)

    # Spearman correlation: layers_active vs centroid_distance
    # Expect negative: more layers -> closer to target centroid
    valid = ~np.isnan(layers) & ~np.isnan(distances)
    if valid.sum() < 10:
        log.warning("Too few valid entries for correlation")
        return pd.DataFrame()

    rho, pval = spearmanr(layers[valid], distances[valid])
    log.info("  Spearman(layers_active, centroid_distance) = %.4f (p=%.2e)", rho, pval)

    # Also correlate layers_active with n_active_sources (sanity check)
    n_active = meta["n_active_sources"].values.astype(float)
    rho2, pval2 = spearmanr(layers[valid], n_active[valid])
    log.info("  Spearman(layers_active, n_active_4sources) = %.4f (p=%.2e)", rho2, pval2)

    return pd.DataFrame(
        [
            {
                "comparison": "layers_active_vs_centroid_distance",
                "spearman_rho": rho,
                "pval": pval,
                "n": int(valid.sum()),
            },
            {
                "comparison": "layers_active_vs_n_active_sources",
                "spearman_rho": rho2,
                "pval": pval2,
                "n": int(valid.sum()),
            },
        ]
    )


# ---------------------------------------------------------------------------
# Downstream analysis
# ---------------------------------------------------------------------------


def run_clustering(
    embeddings: np.ndarray, meta: pd.DataFrame, source_names: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """K-means clustering with silhouette evaluation."""
    log.info("=" * 40)
    log.info("Running K-means clustering (k=%d-%d)", K_RANGE.start, K_RANGE.stop - 1)

    best_k = 6
    best_sil = -1.0
    results = []

    for k in K_RANGE:
        km = KMeans(n_clusters=k, random_state=SEED, n_init=10, max_iter=300)
        labels = km.fit_predict(embeddings)
        sil = silhouette_score(embeddings, labels)
        results.append(
            {"method": "Contrastive", "k": k, "silhouette": sil, "inertia": km.inertia_}
        )
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
            "  Cluster %d: %d genes, %d drug targets, %d pos controls, mean_sources=%.2f",
            c,
            n_genes,
            n_drug,
            n_pc,
            mean_sources,
        )

    return cluster_df, pd.DataFrame(results)


def run_umap(embeddings: np.ndarray, meta: pd.DataFrame) -> pd.DataFrame:
    """UMAP dimensionality reduction."""
    import umap

    log.info("Computing UMAP (n=%d genes)", len(embeddings))
    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=30,
        min_dist=0.3,
        metric="euclidean",
        random_state=SEED,
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


def target_centroid_ranking(
    embeddings: np.ndarray, meta: pd.DataFrame
) -> pd.DataFrame:
    """Rank all genes by distance to drug target centroid."""
    is_target = meta["is_drug_target"].values | meta["is_positive_control"].values
    n_targets = is_target.sum()
    log.info("Computing target centroid from %d known target/control genes", n_targets)

    if n_targets == 0:
        log.warning("No drug targets found in atlas; skipping centroid ranking")
        return pd.DataFrame()

    centroid = embeddings[is_target].mean(axis=0, keepdims=True)
    distances = cdist(embeddings, centroid, metric="euclidean").flatten()

    ranking_df = meta[
        [
            "ensembl_id",
            "human_symbol",
            "gene_biotype",
            "n_active_sources",
            "is_drug_target",
            "is_positive_control",
        ]
    ].copy()
    ranking_df["centroid_distance"] = distances
    ranking_df["rank"] = (
        ranking_df["centroid_distance"].rank(method="min").astype(int)
    )
    ranking_df = ranking_df.sort_values("centroid_distance")

    # Log top 20 novel genes
    novel = ranking_df[
        ~ranking_df["is_drug_target"] & ~ranking_df["is_positive_control"]
    ]
    log.info("Top 20 novel genes nearest drug target centroid:")
    for _, row in novel.head(20).iterrows():
        log.info(
            "  Rank %4d: %-12s dist=%.4f  sources=%d",
            row["rank"],
            row["human_symbol"],
            row["centroid_distance"],
            row["n_active_sources"],
        )

    return ranking_df


def source_attention_analysis(
    attn_weights: np.ndarray,
    cluster_labels: np.ndarray,
    source_names: list[str],
) -> pd.DataFrame:
    """Analyze which sources dominate per cluster."""
    n_clusters = len(set(cluster_labels))

    rows = []
    for c in range(n_clusters):
        idx = cluster_labels == c
        mean_attn = attn_weights[idx].mean(axis=0)
        row = {"cluster": c, "n_genes": int(idx.sum())}
        for i, name in enumerate(source_names):
            row[f"attn_{name}"] = float(mean_attn[i])
        row["dominant_source"] = source_names[int(np.argmax(mean_attn))]
        rows.append(row)

    summary = pd.DataFrame(rows)
    log.info("Source attention by cluster:")
    for _, row in summary.iterrows():
        log.info(
            "  Cluster %d (%d genes): dominant=%s",
            row["cluster"],
            row["n_genes"],
            row["dominant_source"],
        )
    return summary


def run_go_enrichment(
    cluster_df: pd.DataFrame, gene_sets: dict[str, set[str]]
) -> pd.DataFrame:
    """Fisher's exact test Hallmark enrichment for each cluster."""
    log.info("=" * 40)
    log.info("Running GO/Hallmark enrichment per cluster (Fisher's exact test)")

    background = set(cluster_df["human_symbol"].dropna().unique())
    all_rows = []

    for cluster_id in sorted(cluster_df["cluster"].unique()):
        cluster_genes = set(
            cluster_df.loc[
                cluster_df["cluster"] == cluster_id, "human_symbol"
            ]
            .dropna()
            .unique()
        )
        enrich = fisher_enrichment(cluster_genes, background, gene_sets, min_overlap=3)
        if enrich.empty:
            log.info("  Cluster %d: no enrichments (min_overlap=3)", cluster_id)
            continue

        enrich.insert(0, "cluster", cluster_id)
        sig = enrich[enrich["padj"] < 0.05]
        log.info(
            "  Cluster %d: %d sig gene sets (padj<0.05) of %d tested",
            cluster_id,
            len(sig),
            len(enrich),
        )
        if len(sig) > 0:
            for _, row in sig.head(3).iterrows():
                log.info(
                    "    %s: OR=%.1f, padj=%.2e, overlap=%d",
                    row["gene_set"],
                    row["odds_ratio"],
                    row["padj"],
                    row["overlap"],
                )
        all_rows.append(enrich)

    if all_rows:
        return pd.concat(all_rows, ignore_index=True)
    return pd.DataFrame()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    t0 = time.time()
    log.info("=" * 70)
    log.info("Script 171v2: Contrastive Multi-Source Gene Embeddings (4 real sources)")
    log.info("=" * 70)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    log.info("Output directory: %s", RESULTS_DIR)

    # 1. Load and prepare data (4 independent sources only)
    source_arrays, source_masks, source_dims, source_names, meta = (
        load_and_prepare_data()
    )
    log.info("Sources: %s", source_names)
    log.info("Source dims: %s", source_dims)

    # 2. PCA baseline (before training, for fair comparison)
    pca_emb, pca_best_sil, pca_eval = run_pca_baseline(source_arrays, source_masks)

    # 3. Train contrastive model
    model, embeddings, attn_weights, training_metrics = train_model(
        source_arrays, source_masks, source_dims, meta
    )

    # 4. UMAP
    umap_df = run_umap(embeddings, meta)

    # 5. Clustering (contrastive embeddings)
    cluster_df, cluster_eval = run_clustering(embeddings, meta, source_names)

    # 6. Combine baseline comparison
    baseline_comparison = pd.concat([pca_eval, cluster_eval], ignore_index=True)
    pca_best = pca_eval.loc[pca_eval["silhouette"].idxmax()]
    contr_best = cluster_eval.loc[cluster_eval["silhouette"].idxmax()]
    log.info("=" * 40)
    log.info("COMPARISON: PCA best silhouette=%.4f (k=%d) vs Contrastive=%.4f (k=%d)",
             pca_best["silhouette"], int(pca_best["k"]),
             contr_best["silhouette"], int(contr_best["k"]))

    # 7. Target centroid ranking
    ranking_df = target_centroid_ranking(embeddings, meta)

    # 8. layers_active baseline
    layers_corr = run_layers_active_baseline(embeddings, meta)

    # 9. Source attention analysis
    attn_summary = source_attention_analysis(
        attn_weights, cluster_df["cluster"].values, source_names
    )

    # 10. GO/Hallmark enrichment per cluster
    gene_sets = {}
    if GMT_PATH.exists():
        gene_sets = load_hallmark_gmt(GMT_PATH)
    else:
        log.warning("GMT file not found at %s; skipping GO enrichment", GMT_PATH)

    go_enrichment_df = pd.DataFrame()
    if gene_sets:
        go_enrichment_df = run_go_enrichment(cluster_df, gene_sets)

    # -------------------------------------------------------------------
    # Save all outputs
    # -------------------------------------------------------------------
    log.info("=" * 40)
    log.info("Saving outputs to %s", RESULTS_DIR)

    # Embeddings
    emb_df = pd.DataFrame(
        embeddings,
        columns=[f"emb_{i}" for i in range(embeddings.shape[1])],
    )
    emb_df.insert(0, "ensembl_id", meta["ensembl_id"].values)
    emb_df.insert(1, "human_symbol", meta["human_symbol"].values)
    emb_df.to_csv(RESULTS_DIR / "gene_embeddings_32dim.csv", index=False)
    log.info("  gene_embeddings_32dim.csv: %s", emb_df.shape)

    # Clusters
    cluster_df.to_csv(RESULTS_DIR / "gene_clusters.csv", index=False)
    log.info("  gene_clusters.csv: %s", cluster_df.shape)

    # Cluster evaluation + PCA baseline
    baseline_comparison.to_csv(RESULTS_DIR / "baseline_comparison.csv", index=False)
    log.info("  baseline_comparison.csv: %d rows", len(baseline_comparison))

    # Target ranking
    if not ranking_df.empty:
        ranking_df.to_csv(RESULTS_DIR / "target_centroid_ranking.csv", index=False)
        log.info("  target_centroid_ranking.csv: %s", ranking_df.shape)

    # Attention weights (per gene)
    attn_df = pd.DataFrame(
        attn_weights, columns=[f"attn_{n}" for n in source_names]
    )
    attn_df.insert(0, "ensembl_id", meta["ensembl_id"].values)
    attn_df.insert(1, "human_symbol", meta["human_symbol"].values)
    attn_df.to_csv(RESULTS_DIR / "source_attention_weights.csv", index=False)
    log.info("  source_attention_weights.csv: %s", attn_df.shape)

    # Attention summary per cluster
    attn_summary.to_csv(RESULTS_DIR / "source_attention_by_cluster.csv", index=False)

    # Training metrics
    training_metrics.to_csv(RESULTS_DIR / "training_metrics.csv", index=False)
    log.info("  training_metrics.csv: %d epochs", len(training_metrics))

    # UMAP coordinates
    umap_df.to_csv(RESULTS_DIR / "umap_coordinates.csv", index=False)
    log.info("  umap_coordinates.csv: %s", umap_df.shape)

    # layers_active correlation
    if not layers_corr.empty:
        layers_corr.to_csv(
            RESULTS_DIR / "layers_active_correlation.csv", index=False
        )
        log.info("  layers_active_correlation.csv: %d rows", len(layers_corr))

    # GO enrichment
    if not go_enrichment_df.empty:
        go_enrichment_df.to_csv(
            RESULTS_DIR / "cluster_go_enrichment.csv", index=False
        )
        log.info(
            "  cluster_go_enrichment.csv: %d rows (%d sig at padj<0.05)",
            len(go_enrichment_df),
            (go_enrichment_df["padj"] < 0.05).sum(),
        )

    # Save model
    torch.save(model.state_dict(), RESULTS_DIR / "contrastive_model.pt")
    log.info("  contrastive_model.pt saved")

    # -------------------------------------------------------------------
    # Final summary
    # -------------------------------------------------------------------
    elapsed = time.time() - t0
    log.info("=" * 70)
    log.info("Script 171v2 COMPLETE in %.1f min", elapsed / 60)
    log.info("  Sources: %d (S1 human bulk, S2 mouse bulk, S3 genetic causal, S4 essentiality)", len(source_names))
    log.info("  Embedding dim: %d", EMBED_DIM)
    log.info("  PCA baseline silhouette: %.4f", pca_best["silhouette"])
    log.info("  Contrastive best silhouette: %.4f", contr_best["silhouette"])
    if not go_enrichment_df.empty:
        n_sig = (go_enrichment_df["padj"] < 0.05).sum()
        log.info("  Hallmark enrichments (padj<0.05): %d", n_sig)
    log.info("=" * 70)


if __name__ == "__main__":
    main()
