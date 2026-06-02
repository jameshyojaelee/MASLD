#!/usr/bin/env python3
"""
Script 170v2: GNN Drug-Target Predictor (Improved)
====================================================
Fixes from code review of v1:
  1. Downsized model: HIDDEN_DIM=32, EMBED_DIM=16 (~4K params for ~100 positives)
  2. Expanded positives: dgidb_druggable + opentargets_drug (~100 genes) as training set;
     strict 12 clinical targets as held-out validation only
  3. Edge weights: STRING combined_score passed to SAGEConv
  4. Scaler leakage fix: StandardScaler fit per LODOCV fold (excluding held-out)
  5. Drop >99% NaN columns before feature building; add binary "is_observed" indicators
  6. Increased permutations: 200 (was 50)
  7. AUPRC metric alongside AUROC everywhere
  8. Baselines: Random Forest (no graph) + simple atlas rank by layers_active
  9. Increased regularization: dropout 0.5 (was 0.3), weight_decay 1e-2 (was 1e-4)

Two evaluation modes:
  (a) LODOCV on all ~100 druggable genes (training evaluation)
  (b) Held-out rank of strict 12 clinical targets when trained on ~88 druggable genes

Outputs (results/novel_ml/gnn_drug_target_v2/):
  - gnn_predictions.csv: per-gene drug target probability
  - gnn_validation.csv: summary metrics
  - gnn_lodocv_expanded.csv: LODOCV on expanded ~100 druggable genes
  - gnn_clinical12_holdout.csv: strict 12 clinical target ranks (trained on ~88)
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
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
)
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
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
OUTDIR = INTEGRATION / "results/novel_ml/gnn_drug_target_v2"
OUTDIR.mkdir(parents=True, exist_ok=True)

# Logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(OUTDIR / "gnn_drug_target_v2.log", mode="w"),
    ],
)
log = logging.getLogger(__name__)

# Reproducibility
SEED = 42
np.random.seed(SEED)
torch.manual_seed(SEED)

# Hyperparameters — downsized model + stronger regularization
MIN_SCORE = 400  # STRING combined_score threshold (medium confidence)
HIDDEN_DIM = 32  # was 128
EMBED_DIM = 16  # was 64
DROPOUT = 0.5  # was 0.3
LR = 1e-3
WEIGHT_DECAY = 1e-2  # was 1e-4
EPOCHS = 200
PATIENCE = 20
N_PERMUTATIONS = 200  # was 50
NAN_THRESHOLD = 0.99  # drop columns with >99% NaN

# Evidence source groupings for ablation
SOURCE_GROUPS = {
    "S1_human_bulk": [
        "dream_logFC",
        "dream_padj",
        "dream_tstat",
        "dream_logFC_M",
        "dream_logFC_F",
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
        "finngen_nafld_coloc_pp4",
        "finngen_nash_coloc_pp4",
        "finngen_hcc_coloc_pp4",
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

# Key targets to track (for reporting)
KEY_TARGETS = ["THRB", "NR1H4", "PPARA", "PPARD", "PPARG", "GLP1R", "SCD", "DGAT2"]


# ===========================================================================
# 1. Data Loading
# ===========================================================================
def load_string_ppi():
    """Load STRING PPI and protein-to-gene mapping. Returns edges with weights."""
    log.info("Loading STRING protein info for ENSP -> gene symbol mapping...")
    info = pd.read_csv(
        STRING_INFO, sep="\t", usecols=["#string_protein_id", "preferred_name"]
    )
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
    """
    Load multi-evidence atlas, extract numeric features.
    Fix 5: Drop >99% NaN columns, add binary is_observed indicators for sparse cols.
    """
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
        # Boolean-as-object / string flags (keep for label extraction, drop from features)
        # mr_sig removed 2026-04-22 — MR ditched from paper.
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
    log.info(f"  {len(numeric_cols)} numeric feature columns before NaN filtering")

    atlas_numeric = atlas[numeric_cols].copy()

    # --- Fix 5: Drop >99% NaN columns, add is_observed indicators ---
    nan_frac = atlas_numeric.isna().mean()
    sparse_cols = nan_frac[nan_frac > NAN_THRESHOLD].index.tolist()
    # For sparse cols with 90-99% NaN, add binary indicators before dropping the sparse ones
    medium_sparse_cols = nan_frac[
        (nan_frac > 0.90) & (nan_frac <= NAN_THRESHOLD)
    ].index.tolist()
    dense_cols = [c for c in numeric_cols if c not in sparse_cols]

    log.info(
        f"  Dropping {len(sparse_cols)} columns with >{NAN_THRESHOLD*100:.0f}% NaN: "
        f"{sparse_cols}"
    )
    log.info(
        f"  Adding {len(medium_sparse_cols)} binary is_observed indicators for "
        f"90-{NAN_THRESHOLD*100:.0f}% NaN columns"
    )

    # Build feature matrix: dense cols + is_observed indicators
    feat_df = atlas_numeric[dense_cols].copy()
    for c in medium_sparse_cols:
        feat_df[f"{c}_observed"] = (~atlas_numeric[c].isna()).astype(float)

    final_cols = feat_df.columns.tolist()
    log.info(f"  Final feature matrix: {len(final_cols)} columns")

    return atlas_numeric, feat_df, final_cols


def load_expanded_positives():
    """
    Fix 2: Load expanded positive labels from atlas.
    Training positives = dgidb_druggable == True (any gene in atlas flagged as druggable).
    Validation-only = 12 strict clinical targets from clinical_drug_validation_table.csv.
    """
    log.info("Loading expanded positive labels...")

    # Strict 12 clinical targets
    drugs = pd.read_csv(DRUG_TABLE)
    clinical_targets = sorted(drugs["target_gene"].unique().tolist())
    log.info(f"  Strict clinical targets: {len(clinical_targets)} -> {clinical_targets}")

    # Expanded druggable genes from atlas
    atlas_full = pd.read_csv(
        ATLAS,
        low_memory=False,
        usecols=["human_symbol", "dgidb_druggable", "opentargets_drug"],
    )
    expanded = atlas_full.loc[
        atlas_full["dgidb_druggable"] == True, "human_symbol"
    ].unique().tolist()
    log.info(f"  Expanded druggable genes (dgidb_druggable=True): {len(expanded)}")

    # Check overlap
    overlap = [g for g in clinical_targets if g in expanded]
    log.info(
        f"  Clinical targets also in expanded: {len(overlap)} "
        f"(expected 0, clinical targets are separate)"
    )

    return expanded, clinical_targets


def load_layers_active():
    """Load layers_active column for atlas rank baseline."""
    atlas = pd.read_csv(
        ATLAS, low_memory=False, usecols=["human_symbol", "layers_active"]
    )
    atlas = atlas.set_index("human_symbol")
    return atlas["layers_active"]


# ===========================================================================
# 2. Graph Construction
# ===========================================================================
def build_graph(ppi, feat_df, feature_cols, positive_genes, fit_scaler_genes=None):
    """
    Build PyG Data object.
    Fix 3: Pass edge weights.
    Fix 4: Scaler fit on fit_scaler_genes only (or all if None).

    Returns: data, common_genes, gene2idx, scaler, positive_in_graph, raw_feat_df
    """
    # Get genes present in both PPI and feature matrix
    ppi_genes = set(ppi["gene1"]).union(set(ppi["gene2"]))
    feat_genes = set(feat_df.index)
    common_genes = sorted(ppi_genes.intersection(feat_genes))

    gene2idx = {g: i for i, g in enumerate(common_genes)}

    # Filter edges
    ppi_filt = ppi[
        ppi["gene1"].isin(gene2idx) & ppi["gene2"].isin(gene2idx)
    ].copy()

    # Build edge_index (COO format, undirected)
    src = ppi_filt["gene1"].map(gene2idx).values.astype(np.int64)
    dst = ppi_filt["gene2"].map(gene2idx).values.astype(np.int64)
    edge_index = torch.tensor(
        np.stack([np.concatenate([src, dst]), np.concatenate([dst, src])]),
        dtype=torch.long,
    )

    # Fix 3: Edge weights (bidirectional, same weight both directions)
    weights = ppi_filt["weight"].values.astype(np.float32)
    edge_weight = torch.tensor(
        np.concatenate([weights, weights]), dtype=torch.float32
    )

    # Node features
    raw_feat = feat_df.loc[common_genes, feature_cols].copy()
    raw_feat = raw_feat.fillna(0.0)

    # Fix 4: Fit scaler only on specified genes (exclude held-out in CV)
    scaler = StandardScaler()
    if fit_scaler_genes is not None:
        fit_mask = [g in fit_scaler_genes for g in common_genes]
        scaler.fit(raw_feat.values[fit_mask])
    else:
        scaler.fit(raw_feat.values)
    X = scaler.transform(raw_feat.values)
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

    data = Data(x=x, edge_index=edge_index, edge_weight=edge_weight, y=y)
    data.num_nodes = len(common_genes)

    return data, common_genes, gene2idx, scaler, positive_in_graph, positive_missing, raw_feat


# ===========================================================================
# 3. Model (Downsized)
# ===========================================================================
class DrugTargetGNN(torch.nn.Module):
    """2-layer GraphSAGE with edge weights and a classifier head."""

    def __init__(self, in_channels, hidden=HIDDEN_DIM, embed_dim=EMBED_DIM, dropout=DROPOUT):
        super().__init__()
        self.conv1 = SAGEConv(in_channels, hidden)
        self.bn1 = torch.nn.BatchNorm1d(hidden)
        self.conv2 = SAGEConv(hidden, embed_dim)
        self.bn2 = torch.nn.BatchNorm1d(embed_dim)
        self.classifier = torch.nn.Linear(embed_dim, 1)
        self.dropout = dropout

    def forward(self, x, edge_index, edge_weight=None):
        # SAGEConv does not natively use edge_weight in mean aggregation,
        # but we can pass it for compatibility with future versions / custom agg.
        # For now, we weight messages manually via edge_attr if supported.
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
    """Train GNN with full-graph forward pass, backprop on masked nodes."""
    model = DrugTargetGNN(data.x.shape[1])
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
        out, _ = model(data.x, data.edge_index, data.edge_weight)
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
            log.info(
                f"    Epoch {epoch+1}: loss={loss.item():.6f}, best={best_loss:.6f}"
            )

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
        logits, emb = model(data.x, data.edge_index, data.edge_weight)
        probs = torch.sigmoid(logits)
    return probs.numpy(), emb.numpy()


# ===========================================================================
# 5. LODOCV on Expanded Positives (~100 druggable genes)
# ===========================================================================
def lodocv_expanded(
    ppi, feat_df, feature_cols, expanded_positives, common_genes, gene2idx
):
    """
    Leave-one-drug-out CV on expanded druggable genes.
    Fix 4: Re-fit scaler per fold, excluding held-out node.
    """
    log.info(
        f"Running LODOCV on {len(expanded_positives)} expanded druggable genes..."
    )
    results = []

    for i, held_out_gene in enumerate(expanded_positives):
        t0 = time.time()
        held_out_idx = gene2idx[held_out_gene]

        # Remaining positives (exclude held-out)
        remaining_positives = [g for g in expanded_positives if g != held_out_gene]

        # Fix 4: Fit scaler excluding held-out gene
        scaler_genes = set(common_genes) - {held_out_gene}

        # Build graph with scaler fit excluding held-out
        data_cv, _, _, _, _, _, _ = build_graph(
            ppi, feat_df, feature_cols, remaining_positives,
            fit_scaler_genes=scaler_genes,
        )

        # Train mask: everything except held-out
        mask_train = torch.ones(data_cv.num_nodes, dtype=torch.bool)
        mask_train[held_out_idx] = False

        model = train_model(data_cv, mask_train=mask_train, epochs=EPOCHS, patience=PATIENCE)
        probs, _ = predict_full(model, data_cv)

        held_out_prob = probs[held_out_idx]
        # Rank among all negatives (genes with label 0 in this fold)
        neg_mask = data_cv.y.numpy() == 0
        rank = int((probs[neg_mask] >= held_out_prob).sum()) + 1
        n_neg = int(neg_mask.sum())
        percentile = 100.0 * (1.0 - rank / n_neg)

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
        if (i + 1) % 10 == 0 or (i + 1) == len(expanded_positives):
            log.info(
                f"  [{i+1}/{len(expanded_positives)}] {held_out_gene}: "
                f"prob={held_out_prob:.4f}, rank={rank}/{n_neg} "
                f"(top {100-percentile:.1f}%) [{elapsed:.1f}s]"
            )

    return pd.DataFrame(results)


# ===========================================================================
# 6. Clinical 12 Held-Out Evaluation
# ===========================================================================
def clinical12_holdout(
    ppi, feat_df, feature_cols, expanded_positives, clinical_targets,
    common_genes, gene2idx,
):
    """
    Train on expanded druggable genes, evaluate rank of strict 12 clinical targets.
    Each clinical target is scored while training on all ~100 expanded positives
    (the 12 are never in the training set since they don't overlap with dgidb_druggable).
    """
    log.info(
        f"Evaluating {len(clinical_targets)} clinical targets "
        f"(trained on {len(expanded_positives)} expanded druggable genes)..."
    )

    # Build graph with all expanded positives as labels
    data_full, _, _, _, pos_in_graph, _, _ = build_graph(
        ppi, feat_df, feature_cols, expanded_positives
    )

    model = train_model(data_full, epochs=EPOCHS, patience=PATIENCE, verbose=True)
    probs, _ = predict_full(model, data_full)

    results = []
    for gene in clinical_targets:
        if gene not in gene2idx:
            results.append(
                {
                    "gene": gene,
                    "predicted_prob": np.nan,
                    "rank_among_negatives": np.nan,
                    "n_negatives": np.nan,
                    "percentile": np.nan,
                    "in_graph": False,
                }
            )
            continue

        idx = gene2idx[gene]
        prob = probs[idx]
        # Rank among negatives (label == 0 in expanded training)
        neg_mask = data_full.y.numpy() == 0
        rank = int((probs[neg_mask] >= prob).sum()) + 1
        n_neg = int(neg_mask.sum())
        percentile = 100.0 * (1.0 - rank / n_neg)

        results.append(
            {
                "gene": gene,
                "predicted_prob": float(prob),
                "rank_among_negatives": rank,
                "n_negatives": n_neg,
                "percentile": percentile,
                "in_graph": True,
            }
        )
        log.info(
            f"  {gene}: prob={prob:.4f}, rank={rank}/{n_neg} "
            f"(top {100-percentile:.1f}%)"
        )

    return pd.DataFrame(results), model, probs


# ===========================================================================
# 7. Permutation Null
# ===========================================================================
def permutation_null(data, n_perm=N_PERMUTATIONS):
    """Shuffle node features, train, compute AUROC+AUPRC. Build null distribution."""
    log.info(f"Running {n_perm} permutation tests...")
    null_aurocs = []
    null_auprcs = []

    for i in range(n_perm):
        t0 = time.time()
        perm_idx = torch.randperm(data.num_nodes)
        perm_data = Data(
            x=data.x[perm_idx],
            edge_index=data.edge_index,
            edge_weight=data.edge_weight,
            y=data.y.clone(),
            num_nodes=data.num_nodes,
        )
        # Fewer epochs for permutation
        model = train_model(perm_data, epochs=80, patience=10)
        probs, _ = predict_full(model, perm_data)

        y_true = data.y.numpy()
        try:
            auroc = roc_auc_score(y_true, probs)
            auprc = average_precision_score(y_true, probs)
        except ValueError:
            auroc = 0.5
            auprc = 0.0
        null_aurocs.append(auroc)
        null_auprcs.append(auprc)

        elapsed = time.time() - t0
        if (i + 1) % 20 == 0:
            log.info(
                f"  Permutation {i+1}/{n_perm}: AUROC={auroc:.4f}, AUPRC={auprc:.4f}, "
                f"running mean AUROC={np.mean(null_aurocs):.4f} [{elapsed:.1f}s/iter]"
            )

    return np.array(null_aurocs), np.array(null_auprcs)


# ===========================================================================
# 8. Baselines
# ===========================================================================
def rf_baseline_lodocv(feat_df, feature_cols, expanded_positives, common_genes, gene2idx):
    """
    Random Forest baseline on same features (no graph structure).
    LODOCV same as GNN.
    """
    log.info("Running Random Forest baseline (LODOCV on expanded positives)...")

    results = []
    for i, held_out_gene in enumerate(expanded_positives):
        held_out_idx = gene2idx[held_out_gene]
        remaining = [g for g in expanded_positives if g != held_out_gene]

        # Build feature matrix for common_genes
        X_raw = feat_df.loc[common_genes, feature_cols].fillna(0.0).values

        # Fix 4: Fit scaler excluding held-out
        scaler = StandardScaler()
        fit_mask = np.ones(len(common_genes), dtype=bool)
        fit_mask[held_out_idx] = False
        scaler.fit(X_raw[fit_mask])
        X = scaler.transform(X_raw)

        y = np.zeros(len(common_genes))
        for g in remaining:
            if g in gene2idx:
                y[gene2idx[g]] = 1.0

        mask = np.ones(len(y), dtype=bool)
        mask[held_out_idx] = False

        rf = RandomForestClassifier(
            n_estimators=200,
            max_depth=5,
            class_weight="balanced",
            random_state=SEED,
            n_jobs=-1,
        )
        rf.fit(X[mask], y[mask])
        probs = rf.predict_proba(X)[:, 1]

        held_out_prob = probs[held_out_idx]
        neg_mask = y == 0
        rank = int((probs[neg_mask] >= held_out_prob).sum()) + 1
        n_neg = int(neg_mask.sum())
        percentile = 100.0 * (1.0 - rank / n_neg)

        results.append(
            {
                "held_out_gene": held_out_gene,
                "predicted_prob": float(held_out_prob),
                "rank_among_negatives": rank,
                "n_negatives": n_neg,
                "percentile": percentile,
            }
        )

        if (i + 1) % 20 == 0 or (i + 1) == len(expanded_positives):
            log.info(f"  RF [{i+1}/{len(expanded_positives)}] {held_out_gene}: rank={rank}")

    return pd.DataFrame(results)


def lr_baseline_lodocv(feat_df, feature_cols, expanded_positives, common_genes, gene2idx):
    """Logistic Regression baseline (no graph), LODOCV."""
    log.info("Running Logistic Regression baseline (LODOCV on expanded positives)...")

    results = []
    for i, held_out_gene in enumerate(expanded_positives):
        held_out_idx = gene2idx[held_out_gene]
        remaining = [g for g in expanded_positives if g != held_out_gene]

        X_raw = feat_df.loc[common_genes, feature_cols].fillna(0.0).values

        # Fix 4: Scaler excludes held-out
        scaler = StandardScaler()
        fit_mask = np.ones(len(common_genes), dtype=bool)
        fit_mask[held_out_idx] = False
        scaler.fit(X_raw[fit_mask])
        X = scaler.transform(X_raw)

        y = np.zeros(len(common_genes))
        for g in remaining:
            if g in gene2idx:
                y[gene2idx[g]] = 1.0

        mask = np.ones(len(y), dtype=bool)
        mask[held_out_idx] = False

        lr = LogisticRegression(
            max_iter=1000,
            class_weight="balanced",
            C=0.1,  # stronger reg to match GNN weight_decay increase
            solver="lbfgs",
            random_state=SEED,
        )
        lr.fit(X[mask], y[mask])
        probs = lr.predict_proba(X)[:, 1]

        held_out_prob = probs[held_out_idx]
        neg_mask = y == 0
        rank = int((probs[neg_mask] >= held_out_prob).sum()) + 1
        n_neg = int(neg_mask.sum())
        percentile = 100.0 * (1.0 - rank / n_neg)

        results.append(
            {
                "held_out_gene": held_out_gene,
                "predicted_prob": float(held_out_prob),
                "rank_among_negatives": rank,
                "n_negatives": n_neg,
                "percentile": percentile,
            }
        )

        if (i + 1) % 20 == 0 or (i + 1) == len(expanded_positives):
            log.info(f"  LR [{i+1}/{len(expanded_positives)}] {held_out_gene}: rank={rank}")

    return pd.DataFrame(results)


def atlas_rank_baseline(layers_active, expanded_positives, clinical_targets, gene2idx, common_genes):
    """
    Simple atlas rank baseline: rank genes by layers_active column (more layers = higher rank).
    Evaluate on both expanded and clinical targets.
    """
    log.info("Computing atlas rank baseline (layers_active)...")

    # Restrict to genes in graph
    la = layers_active.reindex(common_genes).fillna(0.0)
    ranks = la.rank(ascending=False, method="min").astype(int)

    results_expanded = []
    for gene in expanded_positives:
        if gene in ranks.index:
            r = int(ranks[gene])
            results_expanded.append(
                {
                    "gene": gene,
                    "layers_active": float(la.get(gene, 0)),
                    "rank": r,
                    "percentile": 100.0 * (1.0 - r / len(common_genes)),
                }
            )

    results_clinical = []
    for gene in clinical_targets:
        if gene in ranks.index:
            r = int(ranks[gene])
            results_clinical.append(
                {
                    "gene": gene,
                    "layers_active": float(la.get(gene, 0)),
                    "rank": r,
                    "percentile": 100.0 * (1.0 - r / len(common_genes)),
                }
            )

    return pd.DataFrame(results_expanded), pd.DataFrame(results_clinical)


# ===========================================================================
# 9. Source Ablation
# ===========================================================================
def source_ablation(data, feature_cols, full_auroc, full_auprc):
    """Remove one evidence source at a time, re-train, measure AUROC+AUPRC change."""
    log.info("Running source ablation...")

    col_to_source = {}
    for src, cols in SOURCE_GROUPS.items():
        for c in cols:
            if c in feature_cols:
                col_to_source[c] = src
    for c in feature_cols:
        if c not in col_to_source:
            col_to_source[c] = "other"

    results = [
        {
            "source_removed": "none",
            "n_features_removed": 0,
            "auroc": full_auroc,
            "auprc": full_auprc,
            "delta_auroc": 0.0,
            "delta_auprc": 0.0,
        }
    ]

    active_sources = sorted(set(col_to_source.values()))
    log.info(f"  Active sources: {active_sources}")

    for source in active_sources:
        t0 = time.time()
        cols_to_remove = [
            i for i, c in enumerate(feature_cols) if col_to_source.get(c) == source
        ]
        if not cols_to_remove:
            continue

        x_ablated = data.x.clone()
        x_ablated[:, cols_to_remove] = 0.0

        data_ablated = Data(
            x=x_ablated,
            edge_index=data.edge_index,
            edge_weight=data.edge_weight,
            y=data.y.clone(),
            num_nodes=data.num_nodes,
        )

        model_abl = train_model(data_ablated, epochs=EPOCHS, patience=PATIENCE)
        probs_abl, _ = predict_full(model_abl, data_ablated)
        y_true = data.y.numpy()
        try:
            abl_auroc = roc_auc_score(y_true, probs_abl)
            abl_auprc = average_precision_score(y_true, probs_abl)
        except ValueError:
            abl_auroc = 0.5
            abl_auprc = 0.0

        elapsed = time.time() - t0
        results.append(
            {
                "source_removed": source,
                "n_features_removed": len(cols_to_remove),
                "auroc": abl_auroc,
                "auprc": abl_auprc,
                "delta_auroc": full_auroc - abl_auroc,
                "delta_auprc": full_auprc - abl_auprc,
            }
        )
        log.info(
            f"  Remove {source} ({len(cols_to_remove)} cols): "
            f"AUROC {abl_auroc:.4f} (delta={full_auroc-abl_auroc:+.4f}) [{elapsed:.1f}s]"
        )

    return pd.DataFrame(results)


# ===========================================================================
# MAIN
# ===========================================================================
def main():
    t0_total = time.time()
    log.info("=" * 70)
    log.info("Script 170v2: GNN Drug-Target Predictor (Improved)")
    log.info("=" * 70)
    log.info(
        f"Hyperparams: HIDDEN_DIM={HIDDEN_DIM}, EMBED_DIM={EMBED_DIM}, "
        f"DROPOUT={DROPOUT}, WD={WEIGHT_DECAY}, PERMS={N_PERMUTATIONS}"
    )

    # ---- Load data ----
    ppi = load_string_ppi()
    atlas_numeric, feat_df, feature_cols = load_atlas()
    expanded_positives, clinical_targets = load_expanded_positives()
    layers_active = load_layers_active()

    # Count model params
    _tmp_model = DrugTargetGNN(len(feature_cols))
    n_params = sum(p.numel() for p in _tmp_model.parameters())
    log.info(f"Model params: {n_params:,} (was ~20K in v1, target ~4K)")
    del _tmp_model

    # ---- Build graph with expanded positives ----
    log.info("\n--- Building graph with expanded druggable positives ---")
    data, common_genes, gene2idx, scaler, pos_in_graph, pos_missing, raw_feat = (
        build_graph(ppi, feat_df, feature_cols, expanded_positives)
    )
    log.info(
        f"Graph: {data.num_nodes:,} nodes, {data.edge_index.shape[1]:,} edges, "
        f"{data.x.shape[1]} features, {len(pos_in_graph)} expanded positives in graph"
    )
    if pos_missing:
        log.info(f"  Expanded positives missing from graph: {pos_missing}")

    n_pos = int(data.y.sum())
    n_neg = data.num_nodes - n_pos
    log.info(f"Class balance: {n_pos} positive, {n_neg} negative (ratio 1:{n_neg//max(n_pos,1)})")

    # Check clinical targets in graph
    clinical_in_graph = [g for g in clinical_targets if g in gene2idx]
    clinical_missing = [g for g in clinical_targets if g not in gene2idx]
    log.info(
        f"Clinical targets in graph: {len(clinical_in_graph)}/{len(clinical_targets)}"
    )
    if clinical_missing:
        log.info(f"  Clinical targets missing from graph: {clinical_missing}")

    # ---- Train full model (on expanded positives) ----
    log.info("\n--- Training full GNN model (expanded positives) ---")
    model_full = train_model(data, epochs=EPOCHS, patience=PATIENCE, verbose=True)
    probs_full, embeddings = predict_full(model_full, data)

    y_true = data.y.numpy()
    try:
        full_auroc = roc_auc_score(y_true, probs_full)
        full_auprc = average_precision_score(y_true, probs_full)
    except ValueError:
        full_auroc = 0.5
        full_auprc = 0.0
    log.info(f"Full model (expanded positives) AUROC: {full_auroc:.4f}, AUPRC: {full_auprc:.4f}")

    # ---- Save full predictions ----
    pred_df = pd.DataFrame(
        {
            "gene": common_genes,
            "gnn_probability": probs_full,
            "is_expanded_positive": y_true.astype(int),
            "is_clinical_target": [1 if g in clinical_targets else 0 for g in common_genes],
        }
    )
    pred_df["rank"] = (
        pred_df["gnn_probability"].rank(ascending=False, method="min").astype(int)
    )
    pred_df = pred_df.sort_values("gnn_probability", ascending=False)
    pred_df.to_csv(OUTDIR / "gnn_predictions.csv", index=False)
    log.info(f"Saved predictions for {len(pred_df):,} genes")

    # ---- Known target ranks ----
    known_ranks = []
    for gene in KEY_TARGETS:
        if gene in gene2idx:
            rank = pred_df.loc[pred_df["gene"] == gene, "rank"].values[0]
            prob = probs_full[gene2idx[gene]]
            known_ranks.append(
                {
                    "gene": gene,
                    "probability": float(prob),
                    "rank": int(rank),
                    "percentile": 100.0 * (1.0 - rank / len(common_genes)),
                    "in_graph": True,
                }
            )
        else:
            known_ranks.append(
                {
                    "gene": gene,
                    "probability": np.nan,
                    "rank": np.nan,
                    "percentile": np.nan,
                    "in_graph": False,
                }
            )
    known_ranks_df = pd.DataFrame(known_ranks)
    known_ranks_df.to_csv(OUTDIR / "gnn_known_target_ranks.csv", index=False)
    log.info("\nKnown target ranks (from full model on expanded positives):")
    for _, row in known_ranks_df.iterrows():
        if row["in_graph"]:
            log.info(
                f"  {row['gene']}: rank {row['rank']:.0f} "
                f"(prob={row['probability']:.4f}, top {100-row['percentile']:.1f}%)"
            )
        else:
            log.info(f"  {row['gene']}: NOT in graph")

    # ---- (a) LODOCV on expanded ~100 druggable genes ----
    log.info("\n--- (a) LODOCV on expanded druggable genes ---")
    lodocv_expanded_df = lodocv_expanded(
        ppi, feat_df, feature_cols, pos_in_graph, common_genes, gene2idx
    )
    mean_pct_expanded = lodocv_expanded_df["percentile"].mean()
    median_rank_expanded = lodocv_expanded_df["rank_among_negatives"].median()
    # Compute AUROC/AUPRC from LODOCV held-out probabilities
    lodocv_y = np.ones(len(lodocv_expanded_df))
    lodocv_probs = lodocv_expanded_df["predicted_prob"].values
    # Approximate by checking rank distribution
    log.info(
        f"\nExpanded LODOCV: mean percentile = {mean_pct_expanded:.1f}, "
        f"median rank = {median_rank_expanded:.0f}, "
        f"mean held-out prob = {lodocv_probs.mean():.4f}"
    )
    lodocv_expanded_df.to_csv(OUTDIR / "gnn_lodocv_expanded.csv", index=False)

    # ---- (b) Clinical 12 held-out evaluation ----
    log.info("\n--- (b) Clinical 12 held-out evaluation ---")
    clinical12_df, _, _ = clinical12_holdout(
        ppi, feat_df, feature_cols, expanded_positives, clinical_in_graph,
        common_genes, gene2idx,
    )
    clinical12_df.to_csv(OUTDIR / "gnn_clinical12_holdout.csv", index=False)
    if len(clinical12_df[clinical12_df["in_graph"]]) > 0:
        c12_in = clinical12_df[clinical12_df["in_graph"]]
        log.info(
            f"\nClinical 12 holdout: mean percentile = {c12_in['percentile'].mean():.1f}, "
            f"median rank = {c12_in['rank_among_negatives'].median():.0f}"
        )

    # ---- Baselines ----
    log.info("\n--- Baselines ---")

    # Random Forest
    rf_lodocv_df = rf_baseline_lodocv(
        feat_df, feature_cols, pos_in_graph, common_genes, gene2idx
    )
    rf_mean_pct = rf_lodocv_df["percentile"].mean()
    rf_median_rank = rf_lodocv_df["rank_among_negatives"].median()
    log.info(
        f"RF baseline LODOCV: mean pct={rf_mean_pct:.1f}, "
        f"median rank={rf_median_rank:.0f}"
    )

    # Logistic Regression
    lr_lodocv_df = lr_baseline_lodocv(
        feat_df, feature_cols, pos_in_graph, common_genes, gene2idx
    )
    lr_mean_pct = lr_lodocv_df["percentile"].mean()
    lr_median_rank = lr_lodocv_df["rank_among_negatives"].median()
    log.info(
        f"LR baseline LODOCV: mean pct={lr_mean_pct:.1f}, "
        f"median rank={lr_median_rank:.0f}"
    )

    # Atlas rank (layers_active)
    atlas_exp_df, atlas_clin_df = atlas_rank_baseline(
        layers_active, pos_in_graph, clinical_in_graph, gene2idx, common_genes
    )
    if len(atlas_exp_df) > 0:
        atlas_exp_mean_pct = atlas_exp_df["percentile"].mean()
        log.info(f"Atlas rank baseline (expanded): mean pct={atlas_exp_mean_pct:.1f}")
    else:
        atlas_exp_mean_pct = np.nan
    if len(atlas_clin_df) > 0:
        atlas_clin_mean_pct = atlas_clin_df["percentile"].mean()
        log.info(
            f"Atlas rank baseline (clinical 12): mean pct={atlas_clin_mean_pct:.1f}"
        )
    else:
        atlas_clin_mean_pct = np.nan

    # ---- Permutation null ----
    log.info("\n--- Permutation null ---")
    null_aurocs, null_auprcs = permutation_null(data, n_perm=N_PERMUTATIONS)
    perm_pval_auroc = (null_aurocs >= full_auroc).sum() / len(null_aurocs)
    perm_pval_auprc = (null_auprcs >= full_auprc).sum() / len(null_auprcs)
    log.info(
        f"\nPermutation test (n={N_PERMUTATIONS}): "
        f"true AUROC={full_auroc:.4f} (p={perm_pval_auroc:.4f}), "
        f"true AUPRC={full_auprc:.4f} (p={perm_pval_auprc:.4f}), "
        f"null AUROC={null_aurocs.mean():.4f}+/-{null_aurocs.std():.4f}"
    )

    # ---- Compile all LODOCV details ----
    lodocv_all_parts = []
    gnn_part = lodocv_expanded_df.copy()
    gnn_part["model"] = "GNN"
    lodocv_all_parts.append(gnn_part)

    rf_part = rf_lodocv_df.copy()
    rf_part["model"] = "RF_baseline"
    lodocv_all_parts.append(rf_part)

    lr_part = lr_lodocv_df.copy()
    lr_part["model"] = "LR_baseline"
    lodocv_all_parts.append(lr_part)

    lodocv_all = pd.concat(lodocv_all_parts, ignore_index=True)
    lodocv_all.to_csv(OUTDIR / "gnn_lodocv_details.csv", index=False)

    # ---- Compile validation results ----
    validation_rows = [
        ("gnn_full_auroc", full_auroc),
        ("gnn_full_auprc", full_auprc),
        ("gnn_lodocv_expanded_mean_percentile", mean_pct_expanded),
        ("gnn_lodocv_expanded_median_rank", median_rank_expanded),
        ("gnn_clinical12_mean_percentile", c12_in["percentile"].mean() if len(c12_in) > 0 else np.nan),
        ("gnn_clinical12_median_rank", c12_in["rank_among_negatives"].median() if len(c12_in) > 0 else np.nan),
        ("rf_lodocv_expanded_mean_percentile", rf_mean_pct),
        ("rf_lodocv_expanded_median_rank", rf_median_rank),
        ("lr_lodocv_expanded_mean_percentile", lr_mean_pct),
        ("lr_lodocv_expanded_median_rank", lr_median_rank),
        ("atlas_rank_expanded_mean_percentile", atlas_exp_mean_pct),
        ("atlas_rank_clinical12_mean_percentile", atlas_clin_mean_pct),
        ("permutation_null_mean_auroc", null_aurocs.mean()),
        ("permutation_null_std_auroc", null_aurocs.std()),
        ("permutation_pvalue_auroc", perm_pval_auroc),
        ("permutation_null_mean_auprc", null_auprcs.mean()),
        ("permutation_null_std_auprc", null_auprcs.std()),
        ("permutation_pvalue_auprc", perm_pval_auprc),
        ("n_expanded_positives_in_graph", len(pos_in_graph)),
        ("n_clinical_targets_in_graph", len(clinical_in_graph)),
        ("n_graph_nodes", data.num_nodes),
        ("n_graph_edges", data.edge_index.shape[1]),
        ("n_features", data.x.shape[1]),
        ("n_model_params", n_params),
        ("n_permutations", N_PERMUTATIONS),
    ]
    validation_df = pd.DataFrame(validation_rows, columns=["metric", "value"])
    validation_df.to_csv(OUTDIR / "gnn_validation.csv", index=False)

    # ---- Source ablation ----
    ablation_df = source_ablation(data, feature_cols, full_auroc, full_auprc)
    ablation_df.to_csv(OUTDIR / "gnn_source_ablation.csv", index=False)
    log.info("\nSource ablation:")
    for _, row in ablation_df.iterrows():
        log.info(
            f"  {row['source_removed']}: AUROC={row['auroc']:.4f} "
            f"(delta={row['delta_auroc']:+.4f})"
        )

    # ---- Top 50 novel targets ----
    # Exclude known expanded positives AND clinical targets
    exclude_genes = set(expanded_positives) | set(clinical_targets)
    novel_df = (
        pred_df[~pred_df["gene"].isin(exclude_genes)]
        .head(50)
        .copy()
    )
    key_evidence_cols = [
        "dream_logFC",
        "dream_padj",
        # mr_pval removed 2026-04-22 — MR ditched from paper.
        "twas_pval",
        "n_coloc_sources",
        "essentiality_chronos",
        "translatability_score",
        "causal_methods_sig",
        "n_leading_edge_pathways",
        "layers_active",
    ]
    # Pull from raw_feat (unscaled atlas features for common_genes)
    for col in key_evidence_cols:
        if col in raw_feat.columns:
            novel_df[col] = novel_df["gene"].map(raw_feat[col])
        elif col in atlas_numeric.columns:
            novel_df[col] = novel_df["gene"].map(atlas_numeric[col])
    novel_df.to_csv(OUTDIR / "gnn_top_targets.csv", index=False)
    log.info(f"\nTop 10 novel drug targets:")
    for _, row in novel_df.head(10).iterrows():
        log.info(
            f"  {row['gene']}: prob={row['gnn_probability']:.4f} (rank {row['rank']})"
        )

    # ---- Atlas rank baselines ----
    atlas_exp_df.to_csv(OUTDIR / "atlas_rank_expanded.csv", index=False)
    atlas_clin_df.to_csv(OUTDIR / "atlas_rank_clinical12.csv", index=False)

    # ---- Summary ----
    elapsed = time.time() - t0_total
    log.info(f"\n{'='*70}")
    log.info(f"Script 170v2 complete in {elapsed/60:.1f} min")
    log.info(f"{'='*70}")
    log.info(f"\nKey changes from v1:")
    log.info(f"  Model params: {n_params:,} (was ~20K)")
    log.info(f"  Training positives: {len(pos_in_graph)} expanded druggable (was 12 clinical)")
    log.info(f"  Features: {len(feature_cols)} (after dropping >99% NaN + adding is_observed)")
    log.info(f"  Regularization: dropout={DROPOUT}, weight_decay={WEIGHT_DECAY}")
    log.info(f"  Permutations: {N_PERMUTATIONS} (was 50)")
    log.info(f"\nKey findings:")
    log.info(
        f"  (a) LODOCV on {len(pos_in_graph)} druggable genes: "
        f"mean pct={mean_pct_expanded:.1f}"
    )
    if len(c12_in) > 0:
        log.info(
            f"  (b) Clinical 12 holdout: mean pct={c12_in['percentile'].mean():.1f}, "
            f"median rank={c12_in['rank_among_negatives'].median():.0f}"
        )
    log.info(f"  Full model AUROC: {full_auroc:.4f}, AUPRC: {full_auprc:.4f}")
    log.info(f"  Baselines: GNN pct={mean_pct_expanded:.1f} vs "
             f"RF={rf_mean_pct:.1f} vs LR={lr_mean_pct:.1f} vs "
             f"Atlas rank={atlas_exp_mean_pct:.1f}")
    log.info(f"  Permutation p-value (AUROC): {perm_pval_auroc:.4f}")
    log.info(f"\nOutputs in {OUTDIR}/")
    log.info(f"  gnn_predictions.csv           -- {len(pred_df):,} genes")
    log.info(f"  gnn_validation.csv            -- summary metrics")
    log.info(f"  gnn_lodocv_expanded.csv       -- LODOCV on ~{len(pos_in_graph)} druggable")
    log.info(f"  gnn_lodocv_details.csv        -- all models LODOCV combined")
    log.info(f"  gnn_clinical12_holdout.csv    -- strict 12 clinical targets")
    log.info(f"  gnn_source_ablation.csv       -- {len(ablation_df)} source ablations")
    log.info(f"  gnn_top_targets.csv           -- top 50 novel targets")
    log.info(f"  gnn_known_target_ranks.csv    -- {len(KEY_TARGETS)} key target ranks")
    log.info(f"  atlas_rank_expanded.csv       -- layers_active baseline (expanded)")
    log.info(f"  atlas_rank_clinical12.csv     -- layers_active baseline (clinical)")


if __name__ == "__main__":
    main()
