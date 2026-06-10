#!/usr/bin/env python3
"""
Script 170: GNN Drug-Target Predictor
======================================
Graph neural network on STRING PPI annotated with multi-evidence features
to predict MASLD drug targets.

Architecture:
  - GraphSAGE (2-layer) on STRING PPI graph (score > 400)
  - Node features from multi-evidence atlas (95 numeric columns, imputed + scaled)
  - Positive labels: 12 unique gene targets from clinical_drug_validation_table.csv
  - Validation: Leave-one-drug-out CV, permutation null, logistic regression baseline,
    source ablation

Outputs (results/novel_ml/gnn_drug_target/):
  - gnn_predictions.csv: per-gene drug target probability
  - gnn_validation.csv: LODOCV + permutation results
  - gnn_source_ablation.csv: per-source contribution
  - gnn_top_targets.csv: top 50 novel targets with evidence breakdown
  - gnn_known_target_ranks.csv: where THRB, NR1H4, PPARA rank

Environment: micromamba activate spatial
"""

import os
import sys
import time
import logging
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch_geometric.nn import SAGEConv
from torch_geometric.data import Data
from torch_geometric.loader import NeighborLoader
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
)
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
INTEGRATION = PROJECT_ROOT / "RNA-seq/Human/Patient_Cohorts/analysis/integration"
STRING_PPI = PROJECT_ROOT / "data/string_ppi/9606.protein.links.v12.0.txt.gz"
STRING_INFO = PROJECT_ROOT / "data/string_ppi/9606.protein.info.v12.0.txt.gz"
ATLAS = PROJECT_ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
DRUG_TABLE = (
    PROJECT_ROOT
    / "RNA-seq/results/drug_repurposing/clinical_drug_validation_table.csv"
)
OUTDIR = INTEGRATION / "results/novel_ml/gnn_drug_target"
OUTDIR.mkdir(parents=True, exist_ok=True)

# Logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(OUTDIR / "gnn_drug_target.log", mode="w"),
    ],
)
log = logging.getLogger(__name__)

# Reproducibility
SEED = 42
np.random.seed(SEED)
torch.manual_seed(SEED)

# Hyperparameters
MIN_SCORE = 400  # STRING combined_score threshold (medium confidence)
HIDDEN_DIM = 128
EMBED_DIM = 64
DROPOUT = 0.3
LR = 1e-3
WEIGHT_DECAY = 1e-4
EPOCHS = 200
PATIENCE = 20  # early stopping patience
N_PERMUTATIONS = 50  # permutation null iterations
BATCH_SIZE = 2048  # mini-batch size for NeighborLoader
NUM_NEIGHBORS = [15, 10]  # neighbors sampled per layer (2 GNN layers)

# Evidence source groupings for ablation
SOURCE_GROUPS = {
    "S1_human_bulk": [
        "bulk_logFC",
        "bulk_padj",
        "bulk_tstat",
        "bulk_logFC_M",
        "bulk_logFC_F",
        "sex_interaction_padj",
        "nafl_vs_nash_logFC",
        "nafl_vs_nash_padj",
        "nafl_vs_nash_tstat",
        "adv_fib_logFC",
        "adv_fib_padj",
        "nafl_vs_ctrl_logFC",
        "nafl_vs_ctrl_padj",
        "nas_ge5_logFC",
        "nas_ge5_padj",
        "extreme_logFC",
        "extreme_padj",
        "steatosis_ordinal_coef",
        "inflammation_ordinal_coef",
        "ballooning_ordinal_coef",
        "cirrhosis_logFC",
        "cirrhosis_padj",
        "f2_inflection_logFC",
        "f2_inflection_padj",
        "nash_vs_ctrl_logFC",
        "nash_vs_ctrl_padj",
        "early_late_nash_logFC",
        "early_late_nash_padj",
        "nash_vs_nafl_fibadj_logFC",
        "nash_vs_nafl_fibadj_padj",
        "fibrosis_ordinal_coef",
    ],
    "S2_mouse_bulk": [
        "mouse_meta_logFC",
        "mouse_meta_padj",
        "n_diets_sig",
    ],
    "S3_genetic_causal": [
        # mr_beta, mr_pval removed 2026-04-22 — MR ditched from paper.
        "twas_z",
        "twas_pval",
        "coloc_pp4",
        "broadaway_coloc_pp4",
        "ast_coloc_pp4",
        "ggt_coloc_pp4",
        "pdff_coloc_pp4",
        "ukbb_alt_coloc_pp4",
        "n_liver_enzyme_coloc",
        "best_liver_enzyme_pp4",
        "bbj_alt_coloc_pp4",
        "bbj_ast_coloc_pp4",
        "bbj_ggt_coloc_pp4",
        "ghouse_hcc_coloc_pp4",
        "decode_nafl_coloc_pp4",
        "decode_hcc_coloc_pp4",
        "n_coloc_sources",
        "n_ancestry_gwas",
        "n_validated_ancestry_gwas",
        "panukbb_afr_alt_coloc_pp4",
        "panukbb_afr_ast_coloc_pp4",
        "panukbb_afr_ggt_coloc_pp4",
        "panukbb_csa_alt_coloc_pp4",
        "panukbb_csa_ast_coloc_pp4",
        "panukbb_csa_ggt_coloc_pp4",
        # Enhanced MR columns removed 2026-04-22 — MR ditched from paper.
        "otters_broadaway_pval",
        "otters_broadaway_z",
        "otters_n_gwas_sig",
        "ctwas_pip",
        "ctwas_pval",
        "coloc_susie_best_pp4",
        "coloc_susie_n_signals",
        "hyprcoloc_posterior",
        "hyprcoloc_n_traits",
        "causal_methods_sig",
        "ieqtl_interaction_pval",
        "zenodo_coloc_cell_types",
        "zenodo_n_traits_coloc",
    ],
    "S4_essentiality": [
        "essentiality_chronos",
    ],
    "S7_single_cell": [
        "sceqtl_coloc_pp4_hep",
        "sceqtl_coloc_best_pp4",
        "sceqtl_n_cell_types",
    ],
    "derived": [
        "translatability_score",
        "n_concordant_diets",
        "n_leading_edge_pathways",
        "layers_active",
        "pleiotropy_n_traits",
        "pleiotropy_n_domains",
        "progression_twas_n_contrasts",
        "progression_drug_reversal_n",
        "progression_n_contrasts_sig",
        "progression_tau",
    ],
}

# Key targets to track
KEY_TARGETS = ["THRB", "NR1H4", "PPARA", "PPARD", "PPARG", "GLP1R", "SCD", "DGAT2"]


# ===========================================================================
# 1. Data Loading
# ===========================================================================
def load_string_ppi():
    """Load STRING PPI and protein-to-gene mapping."""
    log.info("Loading STRING protein info for ENSP -> gene symbol mapping...")
    info = pd.read_csv(STRING_INFO, sep="\t", usecols=["#string_protein_id", "preferred_name"])
    info.columns = ["string_id", "symbol"]
    ensp_to_symbol = dict(zip(info["string_id"], info["symbol"]))
    log.info(f"  Loaded {len(ensp_to_symbol):,} protein -> symbol mappings")

    log.info(f"Loading STRING PPI (filtering to combined_score > {MIN_SCORE})...")
    ppi = pd.read_csv(STRING_PPI, sep=" ")
    n_raw = len(ppi)
    ppi = ppi[ppi["combined_score"] > MIN_SCORE].copy()
    log.info(f"  {n_raw:,} raw edges -> {len(ppi):,} after score > {MIN_SCORE}")

    # Map protein IDs to gene symbols
    ppi["gene1"] = ppi["protein1"].map(ensp_to_symbol)
    ppi["gene2"] = ppi["protein2"].map(ensp_to_symbol)
    ppi = ppi.dropna(subset=["gene1", "gene2"])

    # Remove self-loops
    ppi = ppi[ppi["gene1"] != ppi["gene2"]]
    log.info(f"  {len(ppi):,} edges after mapping and removing self-loops")

    # Normalize scores to [0,1]
    ppi["weight"] = ppi["combined_score"] / 1000.0

    return ppi[["gene1", "gene2", "weight"]]


def load_atlas():
    """Load multi-evidence atlas and extract numeric features."""
    log.info("Loading multi-evidence atlas...")
    atlas = pd.read_csv(ATLAS, low_memory=False)
    log.info(f"  Atlas shape: {atlas.shape}")

    # Use human_symbol as index
    atlas = atlas.set_index("human_symbol")

    # Drop non-numeric / categorical / identifier columns
    drop_cols = [
        # Identifiers
        "ensembl_id",
        "gene_biotype",
        "mouse_ortholog",
        "mouse_ensembl",
        # Tier/category labels
        "human_consensus_tier",
        "mouse_consensus_tier",
        "primary_category",
        "best_mouse_model",
        "attribution_class",
        "sex_class",
        "pleiotropy_class",
        "pleiotropy_domains",
        "ferroptosis_class",
        "zonation_class",
        "causal_robustness",
        "progression_deconv_class",
        "progression_twas_direction",
        "sex_progression_class",
        "progression_gene_class",
        "progression_peak_contrast",
        # Boolean-as-object / string flags (mr_sig removed 2026-04-22 — MR ditched)
        "is_conserved",
        "is_nafld_specific",
        "is_essential",
        "is_onset_specific",
        "is_progression_specific",
        "prog_is_onset_specific",
        "cross_ancestry_coloc_replication",
        "zenodo_nafld_coloc",
        "ieqtl_disease_interaction",
        "dgidb_druggable",
        "opentargets_drug",
        "lincs_reversal",
        # Cell-type / trait labels
        "ieqtl_cell_type",
        "ukbb_alt_coloc_cell_type",
        "sceqtl_coloc_cell_type",
        "hyprcoloc_traits",
        # Free-text
        "top_pathways",
    ]
    for c in drop_cols:
        if c in atlas.columns:
            atlas = atlas.drop(columns=[c])

    # Keep only numeric columns
    numeric_cols = atlas.select_dtypes(include=[np.number]).columns.tolist()
    log.info(f"  {len(numeric_cols)} numeric feature columns retained")

    atlas_numeric = atlas[numeric_cols].copy()
    return atlas_numeric, numeric_cols


def load_drug_targets():
    """Load positive labels from clinical drug validation table."""
    log.info("Loading clinical drug validation table...")
    drugs = pd.read_csv(DRUG_TABLE)
    targets = drugs["target_gene"].unique().tolist()
    log.info(f"  {len(targets)} unique drug target genes: {targets}")
    return targets


# ===========================================================================
# 2. Graph Construction
# ===========================================================================
def build_graph(ppi, atlas, feature_cols, positive_genes):
    """
    Build PyG Data object:
      - Nodes = genes present in BOTH PPI and atlas
      - Edges = STRING PPI (undirected)
      - Features = atlas numeric columns (imputed + scaled)
      - Labels = 1 for positive drug targets, 0 otherwise
    """
    log.info("Building graph...")

    # Get genes present in both PPI and atlas
    ppi_genes = set(ppi["gene1"]).union(set(ppi["gene2"]))
    atlas_genes = set(atlas.index)
    common_genes = sorted(ppi_genes.intersection(atlas_genes))
    log.info(f"  PPI genes: {len(ppi_genes):,}, Atlas genes: {len(atlas_genes):,}")
    log.info(f"  Common genes (graph nodes): {len(common_genes):,}")

    # Gene to index mapping
    gene2idx = {g: i for i, g in enumerate(common_genes)}

    # Filter edges to common genes
    ppi_filt = ppi[ppi["gene1"].isin(gene2idx) & ppi["gene2"].isin(gene2idx)].copy()
    log.info(f"  Edges in subgraph: {len(ppi_filt):,}")

    # Build edge_index (COO format, undirected — add reverse edges)
    src = ppi_filt["gene1"].map(gene2idx).values.astype(np.int64)
    dst = ppi_filt["gene2"].map(gene2idx).values.astype(np.int64)
    # Bidirectional
    edge_index = torch.tensor(
        np.stack([np.concatenate([src, dst]), np.concatenate([dst, src])]),
        dtype=torch.long,
    )
    log.info(f"  Bidirectional edges: {edge_index.shape[1]:,}")

    # Node features
    feat_df = atlas.loc[common_genes, feature_cols].copy()
    # Impute NaN with 0 (missing evidence = no signal)
    feat_df = feat_df.fillna(0.0)
    # Scale features
    scaler = StandardScaler()
    X = scaler.fit_transform(feat_df.values)
    x = torch.tensor(X, dtype=torch.float32)

    # Labels
    y = torch.zeros(len(common_genes), dtype=torch.float32)
    positive_in_graph = []
    positive_missing = []
    for g in positive_genes:
        if g in gene2idx:
            y[gene2idx[g]] = 1.0
            positive_in_graph.append(g)
        else:
            positive_missing.append(g)
    log.info(
        f"  Positive labels in graph: {len(positive_in_graph)} / {len(positive_genes)}"
    )
    if positive_missing:
        log.info(f"  Missing from graph: {positive_missing}")

    data = Data(x=x, edge_index=edge_index, y=y)
    data.num_nodes = len(common_genes)

    return data, common_genes, gene2idx, scaler, positive_in_graph, feat_df


# ===========================================================================
# 3. Model
# ===========================================================================
class DrugTargetGNN(torch.nn.Module):
    """2-layer GraphSAGE with a classifier head."""

    def __init__(self, in_channels, hidden=128, embed_dim=64, dropout=0.3):
        super().__init__()
        self.conv1 = SAGEConv(in_channels, hidden)
        self.bn1 = torch.nn.BatchNorm1d(hidden)
        self.conv2 = SAGEConv(hidden, embed_dim)
        self.bn2 = torch.nn.BatchNorm1d(embed_dim)
        self.classifier = torch.nn.Linear(embed_dim, 1)
        self.dropout = dropout

    def forward(self, x, edge_index):
        x = self.conv1(x, edge_index)
        x = self.bn1(x)
        x = F.relu(x)
        x = F.dropout(x, p=self.dropout, training=self.training)
        x = self.conv2(x, edge_index)
        x = self.bn2(x)
        embeddings = x
        x = self.classifier(x).squeeze(-1)
        return x, embeddings


# ===========================================================================
# 4. Training utilities
# ===========================================================================
def compute_pos_weight(y):
    """Compute positive class weight for BCE loss."""
    n_pos = y.sum().item()
    n_neg = len(y) - n_pos
    if n_pos == 0:
        return torch.tensor(1.0)
    return torch.tensor(n_neg / max(n_pos, 1))


def train_model(data, mask_train=None, epochs=EPOCHS, patience=PATIENCE, verbose=False):
    """
    Train GNN with full-graph forward pass but only backprop on masked nodes.
    Uses early stopping on training loss (no validation split with only 12 positives).
    """
    model = DrugTargetGNN(
        data.x.shape[1], hidden=HIDDEN_DIM, embed_dim=EMBED_DIM, dropout=DROPOUT
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=10
    )

    if mask_train is None:
        mask_train = torch.ones(data.num_nodes, dtype=torch.bool)

    pos_weight = compute_pos_weight(data.y[mask_train])

    best_loss = float("inf")
    patience_counter = 0
    best_state = None

    for epoch in range(epochs):
        model.train()
        optimizer.zero_grad()
        out, _ = model(data.x, data.edge_index)
        loss = F.binary_cross_entropy_with_logits(
            out[mask_train], data.y[mask_train], pos_weight=pos_weight
        )
        loss.backward()
        optimizer.step()
        scheduler.step(loss.item())

        if loss.item() < best_loss - 1e-5:
            best_loss = loss.item()
            patience_counter = 0
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
        else:
            patience_counter += 1

        if verbose and (epoch + 1) % 50 == 0:
            log.info(f"    Epoch {epoch+1}: loss={loss.item():.6f}, best={best_loss:.6f}")

        if patience_counter >= patience:
            if verbose:
                log.info(f"    Early stopping at epoch {epoch+1}")
            break

    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    return model


def predict_full(model, data):
    """Get probabilities for all nodes using full-graph forward pass."""
    model.eval()
    with torch.no_grad():
        logits, emb = model(data.x, data.edge_index)
        probs = torch.sigmoid(logits)
    return probs.numpy(), emb.numpy()


# ===========================================================================
# 5. Leave-One-Drug-Out CV
# ===========================================================================
def leave_one_drug_out_cv(data, positive_genes, gene2idx):
    """
    For each positive gene, hold it out, train on remaining positives + all negatives,
    predict held-out gene probability.
    """
    log.info("Running Leave-One-Drug-Out CV...")
    results = []

    for i, held_out_gene in enumerate(positive_genes):
        t0 = time.time()
        held_out_idx = gene2idx[held_out_gene]

        # Create train mask: everything except the held-out node
        mask_train = torch.ones(data.num_nodes, dtype=torch.bool)
        mask_train[held_out_idx] = False

        # Modify labels: set held-out to 0 for training
        data_cv = Data(
            x=data.x.clone(),
            edge_index=data.edge_index,
            y=data.y.clone(),
            num_nodes=data.num_nodes,
        )
        data_cv.y[held_out_idx] = 0.0

        model = train_model(data_cv, mask_train=mask_train, epochs=EPOCHS, patience=PATIENCE)
        probs, _ = predict_full(model, data_cv)

        held_out_prob = probs[held_out_idx]
        # Rank among all negatives
        neg_mask = data.y.numpy() == 0
        rank = int((probs[neg_mask] >= held_out_prob).sum()) + 1
        n_neg = int(neg_mask.sum())
        percentile = 100 * (1 - rank / n_neg)

        results.append(
            {
                "held_out_gene": held_out_gene,
                "predicted_prob": float(held_out_prob),
                "rank_among_negatives": rank,
                "n_negatives": n_neg,
                "percentile": percentile,
            }
        )
        elapsed = time.time() - t0
        log.info(
            f"  [{i+1}/{len(positive_genes)}] {held_out_gene}: "
            f"prob={held_out_prob:.4f}, rank={rank}/{n_neg} "
            f"(top {100-percentile:.1f}%) [{elapsed:.1f}s]"
        )

    return pd.DataFrame(results)


# ===========================================================================
# 6. Permutation Null
# ===========================================================================
def permutation_null(data, n_perm=N_PERMUTATIONS):
    """Shuffle node features, train, compute AUROC. Build null distribution."""
    log.info(f"Running {n_perm} permutation tests...")
    null_aurocs = []

    for i in range(n_perm):
        t0 = time.time()
        perm_idx = torch.randperm(data.num_nodes)
        perm_data = Data(
            x=data.x[perm_idx],
            edge_index=data.edge_index,
            y=data.y.clone(),
            num_nodes=data.num_nodes,
        )
        # Fewer epochs for permutation (just need rough convergence)
        model = train_model(perm_data, epochs=80, patience=10)
        probs, _ = predict_full(model, perm_data)

        try:
            auroc = roc_auc_score(data.y.numpy(), probs)
        except ValueError:
            auroc = 0.5
        null_aurocs.append(auroc)

        elapsed = time.time() - t0
        if (i + 1) % 10 == 0:
            log.info(
                f"  Permutation {i+1}/{n_perm}: AUROC={auroc:.4f}, "
                f"running mean={np.mean(null_aurocs):.4f} [{elapsed:.1f}s/iter]"
            )

    return np.array(null_aurocs)


# ===========================================================================
# 7. Logistic Regression Baseline
# ===========================================================================
def lr_baseline(data, positive_genes, gene2idx):
    """
    Logistic regression on same features (no graph structure).
    LODOCV same as GNN.
    """
    log.info("Running logistic regression baseline (LODOCV)...")
    X = data.x.numpy()
    y = data.y.numpy()

    results = []
    for held_out_gene in positive_genes:
        held_out_idx = gene2idx[held_out_gene]

        # Train on everything except held-out
        mask = np.ones(len(y), dtype=bool)
        mask[held_out_idx] = False
        y_train = y.copy()
        y_train[held_out_idx] = 0  # mask the held-out label

        lr = LogisticRegression(
            max_iter=1000,
            class_weight="balanced",
            C=1.0,
            solver="lbfgs",
            random_state=SEED,
        )
        lr.fit(X[mask], y_train[mask])
        probs = lr.predict_proba(X)[:, 1]

        held_out_prob = probs[held_out_idx]
        neg_mask = y == 0
        rank = int((probs[neg_mask] >= held_out_prob).sum()) + 1
        n_neg = int(neg_mask.sum())
        percentile = 100 * (1 - rank / n_neg)

        results.append(
            {
                "held_out_gene": held_out_gene,
                "predicted_prob": float(held_out_prob),
                "rank_among_negatives": rank,
                "n_negatives": n_neg,
                "percentile": percentile,
            }
        )

    # Full-data AUROC for comparison
    lr_full = LogisticRegression(
        max_iter=1000, class_weight="balanced", C=1.0, solver="lbfgs", random_state=SEED
    )
    lr_full.fit(X, y)
    probs_full = lr_full.predict_proba(X)[:, 1]
    try:
        full_auroc = roc_auc_score(y, probs_full)
        full_auprc = average_precision_score(y, probs_full)
    except ValueError:
        full_auroc = 0.5
        full_auprc = 0.0

    return pd.DataFrame(results), full_auroc, full_auprc


# ===========================================================================
# 8. Source Ablation
# ===========================================================================
def source_ablation(data, feature_cols, full_auroc):
    """
    Remove one evidence source at a time, re-train, measure AUROC change.
    """
    log.info("Running source ablation...")

    # Build column-to-source mapping
    col_to_source = {}
    for src, cols in SOURCE_GROUPS.items():
        for c in cols:
            if c in feature_cols:
                col_to_source[c] = src

    # Assign unmatched columns to "other"
    for c in feature_cols:
        if c not in col_to_source:
            col_to_source[c] = "other"

    results = [{"source_removed": "none", "n_features_removed": 0, "auroc": full_auroc, "delta_auroc": 0.0}]

    # Get unique sources that have columns in our features
    active_sources = sorted(set(col_to_source.values()))
    log.info(f"  Active sources: {active_sources}")

    for source in active_sources:
        t0 = time.time()
        # Columns to zero out
        cols_to_remove = [i for i, c in enumerate(feature_cols) if col_to_source.get(c) == source]
        if not cols_to_remove:
            continue

        # Create ablated data (zero out features for this source)
        x_ablated = data.x.clone()
        x_ablated[:, cols_to_remove] = 0.0

        data_ablated = Data(
            x=x_ablated,
            edge_index=data.edge_index,
            y=data.y.clone(),
            num_nodes=data.num_nodes,
        )

        model_abl = train_model(data_ablated, epochs=EPOCHS, patience=PATIENCE)
        probs_abl, _ = predict_full(model_abl, data_ablated)
        try:
            abl_auroc = roc_auc_score(data.y.numpy(), probs_abl)
        except ValueError:
            abl_auroc = 0.5

        delta = full_auroc - abl_auroc
        elapsed = time.time() - t0
        results.append(
            {
                "source_removed": source,
                "n_features_removed": len(cols_to_remove),
                "auroc": abl_auroc,
                "delta_auroc": delta,
            }
        )
        log.info(
            f"  Remove {source} ({len(cols_to_remove)} cols): "
            f"AUROC {abl_auroc:.4f} (delta={delta:+.4f}) [{elapsed:.1f}s]"
        )

    return pd.DataFrame(results)


# ===========================================================================
# MAIN
# ===========================================================================
def main():
    t0 = time.time()
    log.info("=" * 70)
    log.info("Script 170: GNN Drug-Target Predictor")
    log.info("=" * 70)

    # ---- Load data ----
    ppi = load_string_ppi()
    atlas, feature_cols = load_atlas()
    positive_genes = load_drug_targets()

    # ---- Build graph ----
    data, common_genes, gene2idx, scaler, pos_in_graph, feat_df = build_graph(
        ppi, atlas, feature_cols, positive_genes
    )
    log.info(
        f"Graph: {data.num_nodes:,} nodes, {data.edge_index.shape[1]:,} edges, "
        f"{data.x.shape[1]} features, {len(pos_in_graph)} positive labels"
    )

    # Class imbalance info
    n_pos = int(data.y.sum())
    n_neg = data.num_nodes - n_pos
    log.info(f"Class balance: {n_pos} positive, {n_neg} negative (ratio 1:{n_neg//max(n_pos,1)})")

    # ---- Train full model ----
    log.info("\n--- Training full GNN model ---")
    model_full = train_model(data, epochs=EPOCHS, patience=PATIENCE, verbose=True)
    probs_full, embeddings = predict_full(model_full, data)
    try:
        full_auroc = roc_auc_score(data.y.numpy(), probs_full)
        full_auprc = average_precision_score(data.y.numpy(), probs_full)
    except ValueError:
        full_auroc = 0.5
        full_auprc = 0.0
    log.info(f"Full model AUROC: {full_auroc:.4f}, AUPRC: {full_auprc:.4f}")

    # ---- Save full predictions ----
    pred_df = pd.DataFrame(
        {
            "gene": common_genes,
            "gnn_probability": probs_full,
            "is_known_target": data.y.numpy().astype(int),
        }
    )
    # Add rank
    pred_df["rank"] = pred_df["gnn_probability"].rank(ascending=False, method="min").astype(int)
    pred_df = pred_df.sort_values("gnn_probability", ascending=False)
    pred_df.to_csv(OUTDIR / "gnn_predictions.csv", index=False)
    log.info(f"Saved predictions for {len(pred_df):,} genes")

    # ---- Known target ranks ----
    known_ranks = []
    for gene in KEY_TARGETS:
        if gene in gene2idx:
            idx = gene2idx[gene]
            rank = pred_df.loc[pred_df["gene"] == gene, "rank"].values[0]
            prob = probs_full[idx]
            known_ranks.append(
                {
                    "gene": gene,
                    "probability": float(prob),
                    "rank": int(rank),
                    "percentile": 100 * (1 - rank / len(common_genes)),
                    "in_graph": True,
                }
            )
        else:
            known_ranks.append(
                {"gene": gene, "probability": np.nan, "rank": np.nan, "percentile": np.nan, "in_graph": False}
            )
    known_ranks_df = pd.DataFrame(known_ranks)
    known_ranks_df.to_csv(OUTDIR / "gnn_known_target_ranks.csv", index=False)
    log.info(f"\nKnown target ranks:")
    for _, row in known_ranks_df.iterrows():
        if row["in_graph"]:
            log.info(
                f"  {row['gene']}: rank {row['rank']:.0f} "
                f"(prob={row['probability']:.4f}, top {100-row['percentile']:.1f}%)"
            )
        else:
            log.info(f"  {row['gene']}: NOT in graph")

    # ---- LODOCV ----
    lodocv_df = leave_one_drug_out_cv(data, pos_in_graph, gene2idx)
    mean_percentile = lodocv_df["percentile"].mean()
    median_rank = lodocv_df["rank_among_negatives"].median()
    log.info(
        f"\nLODOCV summary: mean percentile = {mean_percentile:.1f}, "
        f"median rank = {median_rank:.0f}"
    )

    # ---- LR Baseline ----
    lr_lodocv_df, lr_auroc, lr_auprc = lr_baseline(data, pos_in_graph, gene2idx)
    lr_mean_pct = lr_lodocv_df["percentile"].mean()
    log.info(
        f"\nLR baseline: AUROC={lr_auroc:.4f}, AUPRC={lr_auprc:.4f}, "
        f"LODOCV mean pct={lr_mean_pct:.1f}"
    )

    # ---- Permutation null ----
    null_aurocs = permutation_null(data, n_perm=N_PERMUTATIONS)
    perm_pval = (null_aurocs >= full_auroc).sum() / len(null_aurocs)
    log.info(
        f"\nPermutation test: true AUROC={full_auroc:.4f}, "
        f"null mean={null_aurocs.mean():.4f} +/- {null_aurocs.std():.4f}, "
        f"p={perm_pval:.4f}"
    )

    # ---- Compile validation results ----
    validation = {
        "metric": [
            "gnn_full_auroc",
            "gnn_full_auprc",
            "gnn_lodocv_mean_percentile",
            "gnn_lodocv_median_rank",
            "lr_full_auroc",
            "lr_full_auprc",
            "lr_lodocv_mean_percentile",
            "permutation_null_mean_auroc",
            "permutation_null_std_auroc",
            "permutation_pvalue",
            "n_positive_labels",
            "n_graph_nodes",
            "n_graph_edges",
            "n_features",
        ],
        "value": [
            full_auroc,
            full_auprc,
            mean_percentile,
            median_rank,
            lr_auroc,
            lr_auprc,
            lr_mean_pct,
            null_aurocs.mean(),
            null_aurocs.std(),
            perm_pval,
            n_pos,
            data.num_nodes,
            data.edge_index.shape[1],
            data.x.shape[1],
        ],
    }
    validation_df = pd.DataFrame(validation)

    # Append LODOCV per-gene details
    lodocv_combined = lodocv_df.copy()
    lodocv_combined["model"] = "GNN"
    lr_lodocv_combined = lr_lodocv_df.copy()
    lr_lodocv_combined["model"] = "LR_baseline"
    lodocv_all = pd.concat([lodocv_combined, lr_lodocv_combined], ignore_index=True)

    validation_df.to_csv(OUTDIR / "gnn_validation.csv", index=False)
    lodocv_all.to_csv(OUTDIR / "gnn_lodocv_details.csv", index=False)
    log.info(f"Saved validation results")

    # ---- Source ablation ----
    ablation_df = source_ablation(data, feature_cols, full_auroc)
    ablation_df.to_csv(OUTDIR / "gnn_source_ablation.csv", index=False)
    log.info(f"\nSource ablation:")
    for _, row in ablation_df.iterrows():
        log.info(
            f"  {row['source_removed']}: AUROC={row['auroc']:.4f} "
            f"(delta={row['delta_auroc']:+.4f})"
        )

    # ---- Top 50 novel targets ----
    # Exclude known targets
    novel_df = pred_df[pred_df["is_known_target"] == 0].head(50).copy()
    # Add key atlas features
    key_evidence_cols = [
        "bulk_logFC",
        "bulk_padj",
        # mr_pval removed 2026-04-22 — MR ditched from paper.
        "twas_pval",
        "coloc_pp4",
        "n_coloc_sources",
        "essentiality_chronos",
        "translatability_score",
        "causal_methods_sig",
        "n_leading_edge_pathways",
    ]
    for col in key_evidence_cols:
        if col in feat_df.columns:
            novel_df[col] = novel_df["gene"].map(feat_df[col])
    novel_df.to_csv(OUTDIR / "gnn_top_targets.csv", index=False)
    log.info(f"\nTop 10 novel drug targets:")
    for _, row in novel_df.head(10).iterrows():
        log.info(f"  {row['gene']}: prob={row['gnn_probability']:.4f} (rank {row['rank']})")

    # ---- Summary ----
    elapsed = time.time() - t0
    log.info(f"\n{'='*70}")
    log.info(f"Script 170 complete in {elapsed/60:.1f} min")
    log.info(f"{'='*70}")
    log.info(f"\nKey findings:")
    log.info(
        f"  GNN AUROC: {full_auroc:.4f} vs LR baseline: {lr_auroc:.4f} "
        f"(graph boost: {full_auroc-lr_auroc:+.4f})"
    )
    log.info(f"  Permutation p-value: {perm_pval:.4f}")
    log.info(
        f"  LODOCV mean percentile: GNN={mean_percentile:.1f} vs LR={lr_mean_pct:.1f}"
    )
    log.info(f"\nOutputs in {OUTDIR}/")
    log.info(f"  gnn_predictions.csv         -- {len(pred_df):,} genes")
    log.info(f"  gnn_validation.csv          -- summary metrics")
    log.info(f"  gnn_lodocv_details.csv      -- per-gene LODOCV")
    log.info(f"  gnn_source_ablation.csv     -- {len(ablation_df)} source ablations")
    log.info(f"  gnn_top_targets.csv         -- top 50 novel targets")
    log.info(f"  gnn_known_target_ranks.csv  -- {len(KEY_TARGETS)} key target ranks")


if __name__ == "__main__":
    main()
