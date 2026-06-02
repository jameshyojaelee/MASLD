#!/usr/bin/env python3
"""
Disease Signature Scoring: Bulk-to-Single-Cell Projection.

Projects bulk RNA-seq NAS/fibrosis stage-specific disease signatures onto the
integrated scRNA-seq atlas (1.23M cells), producing continuous severity scores,
cell-type-specific disease programs, and disease pseudotime trajectories.

Inputs:
    - Bulk dream DEG results per NAS transition (k vs k-1) and fibrosis transition
    - Cumulative signatures (k vs baseline)
    - Per-cell-type pseudobulk DE results
    - Integrated scRNA-seq atlas (h5ad)

Outputs (to results_gpu_v2/disease_signatures/):
    - signature_gene_sets.csv, gene_id_mapping.csv
    - Per-cell severity scores added to atlas .obs
    - Cell-type disease program heatmaps
    - Disease pseudotime trajectories
    - Validation statistics and figures

Usage:
    sbatch run_disease_signature_scoring.sh
"""

import os
import sys
import warnings
import logging
import gc
from pathlib import Path

warnings.filterwarnings("ignore")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
SC_DIR = os.path.join(BASE, "Analysis/SingleCell")
RESULTS = os.path.join(SC_DIR, "results_gpu_v2")
H5AD = os.path.join(RESULTS, "integrated_atlas.h5ad")

DISEASE_SIG_DIR = os.path.join(
    BASE,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures",
)
PSEUDOBULK_DIR = os.path.join(RESULTS, "pseudobulk_de")

OUT_DIR = os.path.join(RESULTS, "disease_signatures")
FIG_DIR = os.path.join(OUT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Optional GPU init (not required — all ops are CPU matrix math)
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.join(SC_DIR, "scripts"))
try:
    from gpu_utils import init_gpu
    init_gpu()
    log.info("GPU initialized")
except Exception:
    log.info("No GPU or gpu_utils unavailable — running CPU-only")

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import h5py
from scipy import stats
from scipy.sparse import issparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# ============================================================================
# Helper: fix h5ad categoricals
# ============================================================================
def fix_ordered_attrs(h5ad_path):
    """Add missing 'ordered' attribute to categoricals in /obs and /var."""
    with h5py.File(h5ad_path, "a") as f:
        for group_name in ["obs", "var"]:
            if group_name not in f:
                continue
            grp = f[group_name]
            for key in grp.keys():
                item = grp[key]
                if isinstance(item, h5py.Group) and "categories" in item:
                    if "ordered" not in item.attrs:
                        item.attrs["ordered"] = False
                        log.info("Fixed missing 'ordered' attr: /%s/%s", group_name, key)


# ============================================================================
# Step 1: Gene ID Mapping & Signature Extraction
# ============================================================================
def build_gene_id_map(adata):
    """Build mapping from gene identifiers to atlas var_names.

    Handles three cases:
    1. Versioned Ensembl IDs (bulk): ENSG00000310526.1 → strip → lookup in gene_ids
    2. Unversioned Ensembl IDs (pseudobulk DE): ENSG00000238009 → direct lookup in gene_ids
    3. Gene symbols: direct match to var_names

    Returns id_to_var dict mapping any input ID to the atlas var_name.
    """
    log.info("Building gene ID mapping (bulk/pseudobulk → atlas var_names)...")

    id_to_var = {}
    var_names_set = set(adata.var_names)

    # 1. Direct var_name match (gene symbols or IDs that ARE var_names)
    for var_name in adata.var_names:
        id_to_var[var_name] = var_name
    log.info("  Direct var_name entries: %d", len(id_to_var))

    # 2. Map via gene_ids column (unversioned Ensembl → var_name)
    gene_id_col = None
    for col in ["gene_ids", "gene_id", "ensembl_id"]:
        if col in adata.var.columns:
            gene_id_col = col
            break

    if gene_id_col:
        n_added = 0
        for var_name, gene_id in zip(adata.var_names, adata.var[gene_id_col]):
            gene_id_str = str(gene_id)
            # Add unversioned Ensembl ID
            unversioned = gene_id_str.split(".")[0]
            if unversioned not in id_to_var:
                id_to_var[unversioned] = var_name
                n_added += 1
            # Also add the full versioned ID if present
            if gene_id_str not in id_to_var:
                id_to_var[gene_id_str] = var_name
        log.info("  Ensembl ID entries from '%s': %d new", gene_id_col, n_added)

    log.info("  Total mapping entries: %d", len(id_to_var))
    return id_to_var


def map_bulk_gene(gene_id, id_to_var):
    """Map a gene ID (versioned/unversioned Ensembl or symbol) to atlas var_name."""
    gene_id_str = str(gene_id)
    # Try exact match first (symbol or full ID)
    if gene_id_str in id_to_var:
        return id_to_var[gene_id_str]
    # Try stripping version
    unversioned = gene_id_str.split(".")[0]
    return id_to_var.get(unversioned)


def extract_stage_signatures(df, contrast_col, gene_col="gene", top_n=100,
                             padj_floor=0.2):
    """Extract UP and DOWN gene sets per stage contrast.

    Args:
        df: DataFrame with logFC, padj, gene, and a contrast column.
        contrast_col: Column name identifying the contrast (e.g., 'contrast').
        gene_col: Column with gene IDs.
        top_n: Max genes per direction per stage.
        padj_floor: Max padj to consider (adaptive floor).

    Returns:
        dict of {contrast_name: {"up": [genes], "down": [genes]}}
    """
    signatures = {}
    for contrast, grp in df.groupby(contrast_col):
        # Apply padj floor filter
        grp_filt = grp[grp["padj"] < padj_floor].copy()

        # Separate UP and DOWN
        up = grp_filt[grp_filt["logFC"] > 0].sort_values("padj").head(top_n)
        down = grp_filt[grp_filt["logFC"] < 0].sort_values("padj").head(top_n)

        signatures[contrast] = {
            "up": up[gene_col].tolist(),
            "down": down[gene_col].tolist(),
        }
        log.info(
            "  %s: %d UP, %d DOWN (padj < %.2f)",
            contrast, len(up), len(down), padj_floor,
        )
    return signatures


def extract_celltype_signatures(pseudobulk_dir, top_n=100, padj_floor=0.2):
    """Extract per-cell-type disease signatures from pseudobulk DE results."""
    signatures = {}
    csv_files = sorted(Path(pseudobulk_dir).glob("*_de.csv"))
    log.info("Found %d pseudobulk DE files", len(csv_files))

    for csv_path in csv_files:
        ct_name = csv_path.stem.replace("_de", "")
        df = pd.read_csv(csv_path)
        filt = df[df["padj"] < padj_floor].copy()

        up = filt[filt["logFC"] > 0].sort_values("padj").head(top_n)
        down = filt[filt["logFC"] < 0].sort_values("padj").head(top_n)

        sig_name = f"ct_{ct_name}_MASLD"
        signatures[sig_name] = {
            "up": up["gene"].tolist(),
            "down": down["gene"].tolist(),
        }
        log.info("  %s: %d UP, %d DOWN", sig_name, len(up), len(down))

    return signatures


def map_signatures_to_atlas(signatures, id_to_var):
    """Map all gene IDs in signatures to atlas var_names.

    Returns mapped signatures and mapping stats.
    """
    mapped_sigs = {}
    stats_rows = []

    for sig_name, directions in signatures.items():
        mapped_sigs[sig_name] = {}
        for direction, genes in directions.items():
            mapped = []
            for g in genes:
                atlas_name = map_bulk_gene(g, id_to_var)
                if atlas_name is not None:
                    mapped.append(atlas_name)
            mapped_sigs[sig_name][direction] = mapped

            n_orig = len(genes)
            n_mapped = len(mapped)
            pct = 100 * n_mapped / n_orig if n_orig > 0 else 0
            stats_rows.append({
                "signature": sig_name,
                "direction": direction,
                "n_original": n_orig,
                "n_mapped": n_mapped,
                "pct_mapped": pct,
            })

    stats_df = pd.DataFrame(stats_rows)
    mean_pct = stats_df["pct_mapped"].mean()
    log.info("Overall mapping: %.1f%% mean success rate", mean_pct)
    return mapped_sigs, stats_df


def signatures_to_net(mapped_sigs):
    """Convert mapped signatures to decoupler net format.

    Each gene set becomes a source; genes are targets with weight +1 (UP) or -1 (DOWN).
    Returns a single combined net and also separate UP/DOWN nets.
    """
    rows_combined = []
    rows_up = []
    rows_down = []

    for sig_name, directions in mapped_sigs.items():
        for gene in directions.get("up", []):
            rows_combined.append({"source": sig_name, "target": gene, "weight": 1.0})
            rows_up.append({"source": f"{sig_name}_UP", "target": gene, "weight": 1.0})
        for gene in directions.get("down", []):
            rows_combined.append({"source": sig_name, "target": gene, "weight": -1.0})
            rows_down.append({"source": f"{sig_name}_DOWN", "target": gene, "weight": 1.0})

    net_combined = pd.DataFrame(rows_combined)
    net_up = pd.DataFrame(rows_up)
    net_down = pd.DataFrame(rows_down)

    # Deduplicate: multiple versioned Ensembl IDs can map to same atlas gene
    net_combined = net_combined.drop_duplicates(subset=["source", "target"])
    net_up = net_up.drop_duplicates(subset=["source", "target"])
    net_down = net_down.drop_duplicates(subset=["source", "target"])

    # Merge UP-only and DOWN-only into one net for independent scoring
    net_directional = pd.concat([net_up, net_down], ignore_index=True)
    net_directional = net_directional.drop_duplicates(subset=["source", "target"])

    return net_combined, net_directional


# ============================================================================
# Step 2: Load Atlas
# ============================================================================
def load_atlas():
    """Load the scRNA-seq atlas with categorical fix."""
    log.info("Fixing h5ad ordered attributes...")
    fix_ordered_attrs(H5AD)

    log.info("Loading atlas: %s", H5AD)
    adata = sc.read_h5ad(H5AD)
    log.info("Loaded %d cells x %d genes", adata.n_obs, adata.n_vars)

    # Use .raw if present (stores log-normalized data)
    if adata.raw is not None:
        log.info("Using adata.raw for expression values")
        adata_expr = adata.raw.to_adata()
    else:
        adata_expr = adata.copy()

    # Check normalization: if max > 50, likely raw counts → normalize
    if issparse(adata_expr.X):
        max_val = adata_expr.X.max()
    else:
        max_val = np.max(adata_expr.X)

    if max_val > 50:
        log.info("Detected raw counts (max=%.0f). Normalizing...", max_val)
        sc.pp.normalize_total(adata_expr, target_sum=1e4)
        sc.pp.log1p(adata_expr)
    else:
        log.info("Data appears pre-normalized (max=%.2f)", max_val)

    # Transfer obs metadata from original adata
    for col in adata.obs.columns:
        if col not in adata_expr.obs.columns:
            adata_expr.obs[col] = adata.obs.loc[adata_expr.obs_names, col].values

    # Transfer obsm (UMAP, PCA, etc.)
    for key in adata.obsm:
        if key not in adata_expr.obsm:
            adata_expr.obsm[key] = adata.obsm[key]

    # Build condition columns
    if "condition_harmonized" in adata_expr.obs.columns:
        adata_expr.obs["condition_coarse"] = adata_expr.obs["condition_harmonized"].values
    elif "condition" in adata_expr.obs.columns:
        adata_expr.obs["condition_coarse"] = adata_expr.obs["condition"].values
    else:
        adata_expr.obs["condition_coarse"] = "Unknown"

    # Finer condition labels using updated MASLD nomenclature
    # NAFLD → MASL (Metabolic-Associated Steatotic Liver)
    # NASH  → MASH (Metabolic-Associated Steatohepatitis)
    # MASLD stays as-is (undifferentiated, datasets without NAFL/NASH split)
    cond_fine_map = {
        "Healthy": "Healthy",
        "NAFLD": "MASL",
        "NASH": "MASH",
        "MASLD": "MASLD",
        "Cirrhotic": "Cirrhotic",
    }
    # Use raw 'condition' column (preserves NAFLD/NASH/Cirrhotic) rather than
    # condition_coarse (which is from condition_harmonized, already collapsed to Healthy/MASLD)
    if "condition" in adata_expr.obs.columns:
        raw_cond = adata_expr.obs["condition"].astype(str)
    else:
        raw_cond = adata_expr.obs["condition_coarse"].astype(str)
    adata_expr.obs["condition_fine"] = raw_cond.map(cond_fine_map).fillna(raw_cond)

    log.info(
        "Condition distribution (coarse):\n%s",
        adata_expr.obs["condition_coarse"].value_counts().to_string(),
    )
    log.info(
        "Condition distribution (fine):\n%s",
        adata_expr.obs["condition_fine"].value_counts().to_string(),
    )

    # Free original to save memory
    del adata
    gc.collect()

    return adata_expr


# ============================================================================
# Step 2b: Pseudobulk Aggregation
# ============================================================================
def create_pseudobulk_adata(adata, sample_col="sample", cell_type_col="cell_type",
                            cell_type_filter=None, min_cells=10):
    """Aggregate per-sample mean expression into pseudobulk AnnData.

    Averages log-normalized expression across all cells per sample, producing
    a small AnnData (n_samples x n_genes) suitable for bulk-like scoring.

    Parameters
    ----------
    cell_type_filter : str, optional
        If set, subset to this cell type before aggregating (e.g., "Hepatocytes").
        Removes composition confounding by scoring a single cell type.
    min_cells : int
        Minimum cells per sample to include (default 10; use 30 for cell-type-filtered).
    """
    label = f"{cell_type_filter} only" if cell_type_filter else "all cells"
    log.info("Creating pseudobulk AnnData (%s, min_cells=%d)...", label, min_cells)

    if sample_col not in adata.obs.columns:
        log.error("No '%s' column in obs — cannot create pseudobulk", sample_col)
        return None

    # Optional cell-type filter
    if cell_type_filter is not None:
        ct_mask = adata.obs[cell_type_col] == cell_type_filter
        n_ct = ct_mask.sum()
        log.info("  Filtering to %s: %d / %d cells", cell_type_filter, n_ct, adata.n_obs)
        if n_ct == 0:
            log.error("  No cells of type '%s' found", cell_type_filter)
            return None
        adata = adata[ct_mask]

    samples = adata.obs[sample_col].astype(str)
    unique_samples = sorted(samples.unique())
    log.info("  %d unique samples", len(unique_samples))

    # Aggregate expression per sample (mean of log-normalized values)
    X = adata.X
    pb_rows = []
    metadata_rows = []

    for sid in unique_samples:
        mask = (samples == sid).values
        n_cells = mask.sum()
        if n_cells < min_cells:
            continue

        # Mean expression (dense — manageable for ~275 samples)
        if issparse(X):
            row_mean = np.asarray(X[mask].mean(axis=0)).flatten()
        else:
            row_mean = X[mask].mean(axis=0)
        pb_rows.append(row_mean)

        # Per-sample metadata (take first cell's metadata)
        first_idx = np.where(mask)[0][0]
        meta = {}
        for col in ["condition_fine", "condition_coarse"]:
            if col in adata.obs.columns:
                meta[col] = str(adata.obs.iloc[first_idx][col])
        meta["n_cells"] = n_cells
        meta["sample"] = sid
        metadata_rows.append(meta)

    if not pb_rows:
        log.error("  No samples passed min_cells=%d threshold", min_cells)
        return None

    pb_X = np.vstack(pb_rows)
    pb_meta = pd.DataFrame(metadata_rows)
    pb_meta.index = pb_meta["sample"]

    pb_adata = ad.AnnData(
        X=pb_X,
        obs=pb_meta,
        var=adata.var.copy(),
    )
    log.info("  Pseudobulk: %d samples x %d genes", pb_adata.n_obs, pb_adata.n_vars)

    # Save cell-type proportions (only for all-cell mode, not filtered)
    if cell_type_filter is None and cell_type_col in adata.obs.columns:
        log.info("  Computing cell-type proportions per sample...")
        ct_counts = adata.obs.groupby([sample_col, cell_type_col]).size().unstack(fill_value=0)
        ct_pct = ct_counts.div(ct_counts.sum(axis=1), axis=0) * 100
        # Add condition
        sample_cond = adata.obs.groupby(sample_col)["condition_fine"].first()
        ct_pct["condition"] = sample_cond
        ct_pct.to_csv(os.path.join(OUT_DIR, "celltype_proportions_per_sample.csv"))
        log.info("  Saved celltype_proportions_per_sample.csv")

    return pb_adata


# ============================================================================
# Step 3: Multi-Method Ensemble Scoring
# ============================================================================
def score_aucell(adata, net_directional, mapped_sigs):
    """Score cells using decoupler AUCell (rank-based enrichment).

    AUCell doesn't support signed weights, so we score UP and DOWN sets separately
    (via net_directional which has {sig}_UP and {sig}_DOWN sources), then compute
    net scores per original signature.
    """
    import decoupler as dc

    log.info("Running AUCell on %d cells, %d directional gene sets...",
             adata.n_obs, net_directional["source"].nunique())
    dc.mt.aucell(adata, net=net_directional, verbose=True)

    # Extract results from obsm
    score_key = "aucell_estimate" if "aucell_estimate" in adata.obsm else "score_aucell"
    acts_dir = adata.obsm[score_key]
    if not isinstance(acts_dir, pd.DataFrame):
        acts_dir = pd.DataFrame(acts_dir, index=adata.obs_names)
    log.info("AUCell directional complete: %d cells x %d sets", acts_dir.shape[0], acts_dir.shape[1])

    # Compute net score per original signature: UP_score - DOWN_score
    net_scores = {}
    for sig_name in mapped_sigs:
        up_col = f"{sig_name}_UP"
        dn_col = f"{sig_name}_DOWN"
        up_vals = acts_dir[up_col] if up_col in acts_dir.columns else 0
        dn_vals = acts_dir[dn_col] if dn_col in acts_dir.columns else 0
        net_score = up_vals - dn_vals
        if isinstance(net_score, (int, float)):
            # Both missing — skip
            continue
        net_scores[sig_name] = net_score

    acts = pd.DataFrame(net_scores)
    log.info("AUCell net scores: %d cells x %d gene sets", acts.shape[0], acts.shape[1])

    # Clean up AUCell obsm to free memory
    for key in list(adata.obsm.keys()):
        if "aucell" in key.lower():
            del adata.obsm[key]
    gc.collect()

    return acts


def score_scanpy(adata, mapped_sigs):
    """Score cells using scanpy score_genes (control-subtracted mean expression)."""
    log.info("Running scanpy score_genes for %d gene sets...", len(mapped_sigs))

    results = {}
    var_names_set = set(adata.var_names)

    for sig_name, directions in mapped_sigs.items():
        for direction in ["up", "down"]:
            gene_list = [g for g in directions.get(direction, []) if g in var_names_set]
            col_name = f"{sig_name}_{direction.upper()}"

            if len(gene_list) < 5:
                log.warning("  %s: only %d genes, skipping score_genes", col_name, len(gene_list))
                results[col_name] = np.full(adata.n_obs, np.nan)
                continue

            ctrl = min(100, len(gene_list))
            try:
                sc.tl.score_genes(adata, gene_list, score_name=col_name, ctrl_size=ctrl)
                results[col_name] = adata.obs[col_name].values.copy()
                # Clean up obs to avoid clutter
                del adata.obs[col_name]
            except Exception as e:
                log.warning("  %s failed: %s", col_name, e)
                results[col_name] = np.full(adata.n_obs, np.nan)

    results_df = pd.DataFrame(results, index=adata.obs_names)
    log.info("scanpy score_genes complete: %d cells x %d scores", *results_df.shape)
    return results_df


def score_zscore(adata, net):
    """Score cells using decoupler z-score enrichment."""
    import decoupler as dc

    log.info("Running z-score on %d cells, %d gene sets...", adata.n_obs, net["source"].nunique())
    dc.mt.zscore(adata, net=net, verbose=True)

    score_key = "zscore_estimate" if "zscore_estimate" in adata.obsm else "score_zscore"
    acts = adata.obsm[score_key]
    if not isinstance(acts, pd.DataFrame):
        acts = pd.DataFrame(acts, index=adata.obs_names)
    log.info("z-score complete: %d cells x %d gene sets", acts.shape[0], acts.shape[1])

    # Clean up z-score obsm to free memory
    for key in list(adata.obsm.keys()):
        if "zscore" in key.lower():
            del adata.obsm[key]
    gc.collect()

    return acts


def zscore_normalize(series):
    """Z-score standardize a series (mean=0, std=1).

    Preserves direction and magnitude of disease signal, unlike rank
    normalization which maps everything to uniform [0,1] and destroys
    the healthy-vs-disease separation.
    """
    valid = series.dropna()
    if len(valid) == 0 or valid.std() == 0:
        return series
    result = pd.Series(np.nan, index=series.index)
    result.loc[valid.index] = (valid - valid.mean()) / valid.std()
    return result


def compute_ensemble_scores(aucell_df, scanpy_df, zscore_df, mapped_sigs):
    """Compute consensus ensemble scores via z-score standardization and averaging.

    For each signature, combine directional scores into net scores,
    z-score standardize each method (mean=0, std=1) to put them on a
    common scale while preserving disease-healthy separation, then average.
    """
    log.info("Computing ensemble consensus scores...")

    # Identify all directional gene set names
    all_sig_names = sorted(mapped_sigs.keys())
    consensus = {}

    for sig_name in all_sig_names:
        method_scores = []

        # --- AUCell: single combined score (net) ---
        if sig_name in aucell_df.columns:
            method_scores.append(("aucell", zscore_normalize(aucell_df[sig_name])))

        # --- scanpy: separate UP and DOWN scores → net ---
        up_col = f"{sig_name}_UP"
        dn_col = f"{sig_name}_DOWN"
        if up_col in scanpy_df.columns or dn_col in scanpy_df.columns:
            up_vals = scanpy_df[up_col] if up_col in scanpy_df.columns else 0
            dn_vals = scanpy_df[dn_col] if dn_col in scanpy_df.columns else 0
            # Handle NaN: if one direction is all NaN, just use the other
            if isinstance(up_vals, pd.Series) and up_vals.isna().all():
                net_scanpy = -dn_vals if isinstance(dn_vals, pd.Series) else pd.Series(0, index=scanpy_df.index)
            elif isinstance(dn_vals, pd.Series) and dn_vals.isna().all():
                net_scanpy = up_vals
            else:
                net_scanpy = up_vals - dn_vals
            if isinstance(net_scanpy, pd.Series):
                method_scores.append(("scanpy", zscore_normalize(net_scanpy)))

        # --- z-score: single combined score (net) ---
        if sig_name in zscore_df.columns:
            method_scores.append(("zscore", zscore_normalize(zscore_df[sig_name])))

        if len(method_scores) == 0:
            log.warning("  %s: no methods produced scores, skipping", sig_name)
            continue

        # Average z-score standardized scores across methods
        score_matrix = pd.concat([s for _, s in method_scores], axis=1)
        consensus[sig_name] = score_matrix.mean(axis=1)

    consensus_df = pd.DataFrame(consensus)
    log.info("Ensemble consensus: %d cells x %d gene sets", *consensus_df.shape)
    return consensus_df


# ============================================================================
# Step 4: Composite Severity Indices
# ============================================================================
def compute_severity_indices(consensus_df, mapped_sigs):
    """Compute NAS and Fibrosis composite severity indices.

    Approach A: Weighted consecutive sum (later stages weighted more).
    Approach B: Cumulative softmax projection (best-matching stage).
    """
    log.info("Computing composite severity indices...")
    n_cells = consensus_df.shape[0]

    results = {}

    # --- Approach A: Weighted Consecutive Sum ---
    for score_type, prefix, max_k in [("NAS", "NAS", 7), ("Fibrosis", "F", 4)]:
        severity = np.zeros(n_cells)
        stages_used = 0

        for k in range(1, max_k + 1):
            if score_type == "NAS":
                sig_name = f"NAS{k}_vs_NAS{k-1}"
            else:
                sig_name = f"F{k}_vs_F{k-1}"

            if sig_name in consensus_df.columns:
                severity += k * consensus_df[sig_name].fillna(0).values
                stages_used += 1

        if stages_used > 0:
            results[f"{score_type}_severity_A"] = severity
            log.info("  Approach A %s: %d stages used", score_type, stages_used)
        else:
            log.warning("  Approach A %s: no consecutive signatures found", score_type)

    # --- Approach B: Cumulative Softmax Projection ---
    for score_type, prefix, max_k in [("NAS", "NAS", 7), ("Fibrosis", "F", 4)]:
        # Collect cumulative signature scores (k vs baseline)
        stage_scores = []
        stage_indices = []

        for k in range(1, max_k + 1):
            if score_type == "NAS":
                sig_name = f"NAS{k}_vs_NAS0"
            else:
                sig_name = f"F{k}_vs_F0"

            if sig_name in consensus_df.columns:
                stage_scores.append(consensus_df[sig_name].fillna(0).values)
                stage_indices.append(k)

        if len(stage_scores) < 2:
            log.warning("  Approach B %s: only %d cumulative stages, skipping", score_type, len(stage_scores))
            continue

        # Softmax weighting across stages
        score_matrix = np.column_stack(stage_scores)  # cells x stages
        # Temperature scaling to control softmax sharpness
        temperature = 1.0
        # Numerical stability: subtract row-max before exp (prevents overflow
        # with z-scored inputs that can be negative or large positive)
        score_matrix_scaled = score_matrix / temperature
        score_matrix_scaled = score_matrix_scaled - score_matrix_scaled.max(axis=1, keepdims=True)
        exp_scores = np.exp(score_matrix_scaled)
        softmax_weights = exp_scores / exp_scores.sum(axis=1, keepdims=True)

        # Weighted average of stage indices
        stage_arr = np.array(stage_indices)
        severity = (softmax_weights * stage_arr).sum(axis=1)
        results[f"{score_type}_severity_B"] = severity
        log.info("  Approach B %s: %d stages, median=%.2f", score_type, len(stage_indices), np.median(severity))

    return pd.DataFrame(results, index=consensus_df.index)


# ============================================================================
# Step 4b: Pseudobulk-Level Scoring
# ============================================================================
def score_pseudobulk(pb_adata, net_combined, net_directional, mapped_sigs,
                     prefix="pseudobulk"):
    """Score pseudobulk samples with 3-method ensemble + severity indices.

    Reuses the same scoring methods as cell-level scoring but operates on
    the pseudobulk AnnData (~275 samples instead of 1.23M cells).

    Parameters
    ----------
    prefix : str
        Filename prefix for outputs (default "pseudobulk"; use "hep_pseudobulk"
        for hepatocyte-only to avoid overwriting all-cell results).
    """
    log.info("=== %s scoring (%d samples) ===", prefix, pb_adata.n_obs)

    # 1. Score with all three methods
    log.info("--- %s AUCell ---", prefix)
    pb_aucell = score_aucell(pb_adata, net_directional, mapped_sigs)

    log.info("--- %s scanpy score_genes ---", prefix)
    pb_scanpy = score_scanpy(pb_adata, mapped_sigs)

    log.info("--- %s z-score ---", prefix)
    pb_zscore = score_zscore(pb_adata, net_combined)

    # 2. Ensemble consensus
    log.info("--- %s ensemble consensus ---", prefix)
    pb_consensus = compute_ensemble_scores(pb_aucell, pb_scanpy, pb_zscore, mapped_sigs)
    pb_consensus.to_csv(os.path.join(OUT_DIR, f"{prefix}_consensus.csv"))

    # 3. Severity indices
    log.info("--- %s severity indices ---", prefix)
    pb_severity = compute_severity_indices(pb_consensus, mapped_sigs)

    # Add condition metadata
    for col in ["condition_fine", "condition_coarse"]:
        if col in pb_adata.obs.columns:
            pb_severity[col] = pb_adata.obs.loc[pb_severity.index, col].values

    pb_severity.to_csv(os.path.join(OUT_DIR, f"{prefix}_severity_scores.csv"))
    log.info("%s severity shape: %s", prefix, pb_severity.shape)

    # 4. Validate: Spearman rho vs condition ordinal
    log.info("--- %s validation ---", prefix)
    ordinal_map = {"Healthy": 0, "MASL": 1, "MASLD": 2, "MASH": 3, "Cirrhotic": 4}
    cond_col = "condition_fine" if "condition_fine" in pb_severity.columns else "condition_coarse"
    pb_severity["ordinal"] = pb_severity[cond_col].map(ordinal_map)

    sev_cols = [c for c in pb_severity.columns if "severity" in c]
    val_rows = []
    for col in sev_cols:
        valid = pb_severity.dropna(subset=["ordinal", col])
        if len(valid) > 10:
            rho, p = stats.spearmanr(valid["ordinal"], valid[col])
            val_rows.append({"score": col, "spearman_rho": rho, "pvalue": p, "n": len(valid)})
            log.info("  %s %s: rho=%.3f, p=%.2e, n=%d", prefix, col, rho, p, len(valid))

            # Per-condition means
            for cond in sorted(ordinal_map.keys()):
                cond_vals = valid[valid[cond_col] == cond][col]
                if len(cond_vals) > 0:
                    log.info("    %s: mean=%.3f (n=%d)", cond, cond_vals.mean(), len(cond_vals))

    val_df = pd.DataFrame(val_rows)
    val_df.to_csv(os.path.join(OUT_DIR, f"validation_{prefix}_severity.csv"), index=False)

    # Clean up ordinal column before return
    pb_severity = pb_severity.drop(columns=["ordinal"], errors="ignore")
    return pb_severity


def plot_pseudobulk_calibration(pb_severity_path, fig_dir):
    """Plot pseudobulk severity by condition with Spearman annotation."""
    log.info("Plotting pseudobulk severity calibration...")
    pb = pd.read_csv(pb_severity_path, index_col=0)

    cond_col = "condition_fine" if "condition_fine" in pb.columns else "condition_coarse"
    ordinal_map = {"Healthy": 0, "MASL": 1, "MASLD": 2, "MASH": 3, "Cirrhotic": 4}
    pb["ordinal"] = pb[cond_col].map(ordinal_map)
    pb = pb.dropna(subset=["ordinal"])

    sev_cols = [c for c in pb.columns if "severity" in c]
    if not sev_cols:
        return

    fig, axes = plt.subplots(1, len(sev_cols), figsize=(5 * len(sev_cols), 4))
    if len(sev_cols) == 1:
        axes = [axes]

    cond_order = ["Healthy", "MASL", "MASLD", "MASH", "Cirrhotic"]
    for ax, col in zip(axes, sev_cols):
        groups, positions, labels = [], [], []
        for i, cond in enumerate(cond_order):
            vals = pb[pb[cond_col] == cond][col].dropna()
            if len(vals) > 0:
                groups.append(vals.values)
                positions.append(i)
                labels.append(f"{cond}\n(n={len(vals)})")

        if groups:
            bp = ax.boxplot(groups, positions=positions, widths=0.4,
                            showfliers=False, patch_artist=True,
                            medianprops=dict(color="black", linewidth=0.8))
            colors = ["#1565C0", "#F48FB1", "#AB47BC", "#C2185B", "#880E4F"]
            for patch, pos in zip(bp["boxes"], positions):
                patch.set_facecolor(colors[pos])
                patch.set_alpha(0.6)
                patch.set_linewidth(0.3)
            for w in bp["whiskers"] + bp["caps"]:
                w.set_linewidth(0.3)

            # Strip overlay
            for pos, vals in zip(positions, groups):
                jitter = np.random.uniform(-0.12, 0.12, len(vals))
                ax.scatter(pos + jitter, vals, s=8, alpha=0.4,
                           color=colors[pos], edgecolors="none", zorder=2)

        # Spearman annotation
        valid = pb.dropna(subset=["ordinal", col])
        if len(valid) > 10:
            rho, p = stats.spearmanr(valid["ordinal"], valid[col])
            ax.set_title(f"{col}\nSpearman ρ = {rho:.3f} (p = {p:.2e})", fontsize=9)
        else:
            ax.set_title(col, fontsize=9)

        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_ylabel("Pseudobulk severity score", fontsize=8)

    plt.tight_layout()
    # Derive output name from input path (e.g., hep_pseudobulk_severity_scores → hep_pseudobulk)
    stem = os.path.basename(pb_severity_path).replace("_severity_scores.csv", "")
    out = os.path.join(fig_dir, f"{stem}_severity_calibration.pdf")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved: %s", out)


# ============================================================================
# Step 5: Cell-Type-Specific Disease Programs
# ============================================================================
def compute_celltype_disease_programs(consensus_df, adata, cell_type_col="cell_type"):
    """Score each cell against cell-type-specific signatures and build activity matrix.

    Returns a cell_type × signature heatmap matrix (mean score per group).
    """
    log.info("Computing cell-type disease program activity matrix...")

    # Identify cell-type-specific signatures (prefix "ct_")
    ct_sigs = [c for c in consensus_df.columns if c.startswith("ct_")]
    if not ct_sigs:
        log.warning("No cell-type signatures found in consensus scores")
        return None

    cell_types = adata.obs[cell_type_col].unique()

    # Build activity matrix: mean score per cell type per signature
    activity = pd.DataFrame(index=sorted(cell_types), columns=sorted(ct_sigs), dtype=float)
    for ct in cell_types:
        mask = adata.obs[cell_type_col] == ct
        ct_idx = adata.obs_names[mask]
        common = ct_idx.intersection(consensus_df.index)
        if len(common) == 0:
            continue
        activity.loc[ct] = consensus_df.loc[common, ct_sigs].mean(axis=0).values

    log.info("Activity matrix: %d cell types x %d signatures", *activity.shape)
    return activity


# ============================================================================
# Step 6: Disease Pseudotime Trajectory
# ============================================================================
def compute_pseudotime(adata, severity_df, cell_type_col="cell_type",
                       target_celltype="Hepatocytes", n_comps=15):
    """Compute disease pseudotime for a cell type using severity-anchored DPT.

    Uses the composite severity score as pseudotime ordering, then identifies
    genes correlated with this disease progression axis.

    Returns adata_sub with pseudotime and correlations.
    """
    log.info("Computing disease pseudotime for %s...", target_celltype)

    mask = adata.obs[cell_type_col] == target_celltype
    if mask.sum() < 500:
        log.warning("Only %d %s cells — too few for pseudotime", mask.sum(), target_celltype)
        return None, None

    adata_sub = adata[mask].copy()
    log.info("Subset: %d %s cells", adata_sub.n_obs, target_celltype)

    # Merge severity scores
    common_idx = adata_sub.obs_names.intersection(severity_df.index)
    for col in severity_df.columns:
        adata_sub.obs[col] = np.nan
        adata_sub.obs.loc[common_idx, col] = severity_df.loc[common_idx, col].values

    # Use Approach B severity as pseudotime (best-matching cumulative stage)
    sev_cols = [c for c in severity_df.columns if "severity" in c]

    # Pick the most appropriate severity score for this cell type
    if target_celltype == "Fibroblasts":
        # Fibroblasts → fibrosis trajectory
        ptime_col = "Fibrosis_severity_B" if "Fibrosis_severity_B" in sev_cols else sev_cols[0]
    else:
        # Default to NAS for hepatocytes, macrophages
        ptime_col = "NAS_severity_B" if "NAS_severity_B" in sev_cols else sev_cols[0]

    log.info("Using %s as pseudotime axis for %s", ptime_col, target_celltype)

    # Rank-normalize severity to [0, 1] as pseudotime
    ptime_raw = adata_sub.obs[ptime_col].values
    valid_mask = ~np.isnan(ptime_raw)
    dpt = np.full(len(ptime_raw), np.nan)
    if valid_mask.sum() > 0:
        dpt[valid_mask] = stats.rankdata(ptime_raw[valid_mask]) / valid_mask.sum()
    adata_sub.obs["dpt_pseudotime"] = dpt
    log.info("Pseudotime: median=%.3f, %d valid cells", np.nanmedian(dpt), valid_mask.sum())

    # Compute PCA + neighbors + UMAP for visualization
    try:
        sc.pp.highly_variable_genes(adata_sub, n_top_genes=2000, flavor="seurat_v3",
                                    subset=True)
        sc.pp.pca(adata_sub, n_comps=n_comps)
        sc.pp.neighbors(adata_sub, n_neighbors=30, use_rep="X_pca")
        sc.tl.umap(adata_sub)
    except Exception as e:
        log.warning("UMAP computation failed for %s: %s", target_celltype, e)

    # Correlate pseudotime with all severity scores
    corr_results = []
    valid_dpt = ~np.isnan(dpt)
    for col in sev_cols:
        vals = adata_sub.obs[col].values
        valid = valid_dpt & ~np.isnan(vals)
        if valid.sum() > 100:
            rho, pval = stats.spearmanr(dpt[valid], vals[valid])
            corr_results.append({"score": col, "spearman_rho": rho, "pvalue": pval, "n": valid.sum()})
            log.info("  Pseudotime vs %s: rho=%.3f, p=%.2e", col, rho, pval)

    corr_df = pd.DataFrame(corr_results) if corr_results else None

    # Identify pseudotime-associated genes (top correlating genes)
    log.info("Finding pseudotime-associated genes...")
    gene_corrs = []
    X = adata_sub.X
    if issparse(X):
        X = X.toarray()
    for i, gene in enumerate(adata_sub.var_names):
        expr = X[:, i]
        valid = valid_dpt & ~np.isnan(expr) & (np.std(expr) > 0)
        if valid.sum() > 100:
            rho, _ = stats.spearmanr(dpt[valid], expr[valid])
            gene_corrs.append({"gene": gene, "dpt_spearman": rho})
    gene_corr_df = pd.DataFrame(gene_corrs).sort_values("dpt_spearman", key=abs, ascending=False)
    log.info("Top pseudotime genes: %s", gene_corr_df.head(10)["gene"].tolist())

    return adata_sub, {"dpt_severity_corr": corr_df, "gene_correlations": gene_corr_df}


# ============================================================================
# Step 7: Validation
# ============================================================================
def validate_scores(adata, severity_df, consensus_df, aucell_df, scanpy_df,
                    zscore_df, mapped_sigs, cell_type_col="cell_type"):
    """Run validation checks: label concordance, cross-method agreement, calibration."""
    log.info("=== Running validation ===")
    validation = {}

    # --- 7a. Coarse label concordance ---
    log.info("7a: Label concordance (fine conditions)...")
    # Use condition_fine for finer ordinal trend testing
    cond_col = "condition_fine" if "condition_fine" in adata.obs.columns else "condition_coarse"
    condition_order = {"Healthy": 0, "MASL": 1, "MASLD": 2, "MASH": 3, "Cirrhotic": 4}
    conditions = adata.obs[cond_col].values
    ordinal = np.array([condition_order.get(str(c), -1) for c in conditions])

    concordance_rows = []
    for col in severity_df.columns:
        vals = severity_df.reindex(adata.obs_names)[col].values
        valid = (ordinal >= 0) & ~np.isnan(vals)
        if valid.sum() < 100:
            continue

        # Kruskal-Wallis
        groups = {}
        for cond, ord_val in condition_order.items():
            mask = (ordinal == ord_val) & ~np.isnan(vals)
            if mask.sum() > 10:
                groups[cond] = vals[mask]
        if len(groups) >= 2:
            kw_stat, kw_p = stats.kruskal(*groups.values())
        else:
            kw_stat, kw_p = np.nan, np.nan

        # Jonckheere-Terpstra (approximation via Spearman on ordinals)
        jt_rho, jt_p = stats.spearmanr(ordinal[valid], vals[valid])

        concordance_rows.append({
            "score": col,
            "kruskal_wallis_stat": kw_stat,
            "kw_pvalue": kw_p,
            "jt_spearman_rho": jt_rho,
            "jt_pvalue": jt_p,
            "n_valid": valid.sum(),
        })
        log.info("  %s: KW p=%.2e, trend rho=%.3f (p=%.2e)", col, kw_p, jt_rho, jt_p)

    validation["concordance"] = pd.DataFrame(concordance_rows)

    # --- 7b. Pseudobulk calibration ---
    log.info("7b: Pseudobulk calibration (per-sample aggregation)...")
    if "sample" in adata.obs.columns:
        samples = adata.obs["sample"].values
        sev_reindexed = severity_df.reindex(adata.obs_names)
        sample_scores = sev_reindexed.copy()
        sample_scores["sample"] = samples
        sample_scores["condition"] = conditions
        pseudobulk_scores = sample_scores.groupby("sample").agg(
            {col: "mean" for col in severity_df.columns}
        )
        # Add condition per sample (use fine labels)
        sample_cond = sample_scores.groupby("sample")["condition"].first()
        pseudobulk_scores["condition"] = sample_cond
        validation["pseudobulk_calibration"] = pseudobulk_scores
        log.info("  Pseudobulk scores: %d samples", len(pseudobulk_scores))
    else:
        log.warning("  No 'sample' column — skipping pseudobulk calibration")

    # --- 7c. Cross-method agreement ---
    log.info("7c: Cross-method agreement...")
    agreement_rows = []
    for sig_name in mapped_sigs:
        methods = {}
        if sig_name in aucell_df.columns:
            methods["aucell"] = aucell_df[sig_name]
        if sig_name in zscore_df.columns:
            methods["zscore"] = zscore_df[sig_name]
        # scanpy: compute net from UP-DOWN
        up_col = f"{sig_name}_UP"
        dn_col = f"{sig_name}_DOWN"
        if up_col in scanpy_df.columns or dn_col in scanpy_df.columns:
            up = scanpy_df[up_col] if up_col in scanpy_df.columns else 0
            dn = scanpy_df[dn_col] if dn_col in scanpy_df.columns else 0
            net_sp = up - dn
            if isinstance(net_sp, pd.Series) and not net_sp.isna().all():
                methods["scanpy"] = net_sp

        if len(methods) < 2:
            continue

        method_names = list(methods.keys())
        for i in range(len(method_names)):
            for j in range(i + 1, len(method_names)):
                m1, m2 = method_names[i], method_names[j]
                v1 = methods[m1].dropna()
                v2 = methods[m2].dropna()
                common = v1.index.intersection(v2.index)
                if len(common) > 100:
                    rho, _ = stats.spearmanr(v1[common], v2[common])
                    agreement_rows.append({
                        "signature": sig_name,
                        "method1": m1,
                        "method2": m2,
                        "spearman_rho": rho,
                        "n": len(common),
                    })

    agreement_df = pd.DataFrame(agreement_rows)
    if len(agreement_df) > 0:
        mean_rho = agreement_df["spearman_rho"].mean()
        low_agreement = agreement_df[agreement_df["spearman_rho"] < 0.5]
        log.info("  Mean cross-method rho=%.3f; %d/%d gene sets with rho<0.5",
                 mean_rho, len(low_agreement["signature"].unique()),
                 len(agreement_df["signature"].unique()))
    validation["cross_method_agreement"] = agreement_df

    return validation


# ============================================================================
# Step 8: Visualization
# ============================================================================
def plot_umap_severity(adata, severity_df, fig_dir):
    """UMAP colored by NAS and fibrosis severity scores."""
    log.info("Plotting UMAP severity overlays...")

    if "X_umap" not in adata.obsm:
        log.warning("No X_umap in adata — skipping UMAP plots")
        return

    sev_cols = [c for c in severity_df.columns if "severity" in c]
    if not sev_cols:
        return

    n_panels = len(sev_cols)
    fig, axes = plt.subplots(1, n_panels, figsize=(6 * n_panels, 5))
    if n_panels == 1:
        axes = [axes]

    umap = adata.obsm["X_umap"]
    for ax, col in zip(axes, sev_cols):
        vals = severity_df.reindex(adata.obs_names)[col].values
        valid = ~np.isnan(vals)
        sc_plot = ax.scatter(
            umap[valid, 0], umap[valid, 1],
            c=vals[valid], cmap="viridis", s=0.3, alpha=0.5, rasterized=True,
        )
        plt.colorbar(sc_plot, ax=ax, shrink=0.7)
        ax.set_title(col.replace("_", " "), fontsize=10)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_xticks([])
        ax.set_yticks([])

    plt.tight_layout()
    out = os.path.join(fig_dir, "umap_severity_scores.pdf")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved: %s", out)


def plot_celltype_violins(adata, severity_df, fig_dir, cell_type_col="cell_type"):
    """Per-cell-type violin plots of severity by condition."""
    log.info("Plotting cell-type violin plots...")

    sev_cols = [c for c in severity_df.columns if "severity" in c]
    if not sev_cols:
        return

    cond_col = "condition_fine" if "condition_fine" in adata.obs.columns else "condition_coarse"
    condition_order = ["Healthy", "MASL", "MASLD", "MASH", "Cirrhotic"]
    cell_types = sorted(adata.obs[cell_type_col].unique())

    for score_col in sev_cols:
        plot_df = pd.DataFrame({
            "cell_type": adata.obs[cell_type_col].values,
            "condition": adata.obs[cond_col].values,
            "score": severity_df.reindex(adata.obs_names)[score_col].values,
        })
        plot_df = plot_df[plot_df["condition"].isin(condition_order)].dropna(subset=["score"])

        n_ct = len(cell_types)
        fig, axes = plt.subplots(1, n_ct, figsize=(3 * n_ct, 4), sharey=True)
        if n_ct == 1:
            axes = [axes]

        for ax, ct in zip(axes, cell_types):
            ct_data = plot_df[plot_df["cell_type"] == ct]
            if len(ct_data) < 10:
                ax.set_title(ct, fontsize=8)
                continue

            positions = []
            data_groups = []
            labels = []
            for i, cond in enumerate(condition_order):
                vals = ct_data[ct_data["condition"] == cond]["score"].values
                if len(vals) > 0:
                    data_groups.append(vals)
                    positions.append(i)
                    labels.append(cond)

            if data_groups:
                parts = ax.violinplot(data_groups, positions=positions, showmedians=True)
                for pc in parts.get("bodies", []):
                    pc.set_alpha(0.7)
                ax.set_xticks(positions)
                ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)

            ax.set_title(ct, fontsize=8)
            if ax == axes[0]:
                ax.set_ylabel(score_col)

        plt.suptitle(score_col.replace("_", " "), fontsize=12)
        plt.tight_layout()
        out = os.path.join(fig_dir, f"violins_{score_col}.pdf")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        log.info("Saved: %s", out)


def plot_celltype_stage_heatmap(consensus_df, adata, fig_dir, cell_type_col="cell_type"):
    """Heatmap: mean consensus score per cell type per stage signature."""
    log.info("Plotting cell-type × stage heatmap...")

    # Filter to stage signatures only (NAS and fibrosis consecutive)
    stage_cols = [c for c in consensus_df.columns
                  if ("NAS" in c and "vs" in c) or ("F" in c and "vs" in c)]
    stage_cols = [c for c in stage_cols if not c.startswith("ct_")]

    if not stage_cols:
        log.warning("No stage signatures for heatmap")
        return

    cell_types = sorted(adata.obs[cell_type_col].unique())
    heatmap = pd.DataFrame(index=cell_types, columns=stage_cols, dtype=float)

    for ct in cell_types:
        mask = adata.obs[cell_type_col] == ct
        ct_idx = adata.obs_names[mask].intersection(consensus_df.index)
        if len(ct_idx) > 0:
            heatmap.loc[ct] = consensus_df.loc[ct_idx, stage_cols].mean(axis=0).values

    fig, ax = plt.subplots(figsize=(max(12, len(stage_cols) * 0.6), max(6, len(cell_types) * 0.4)))
    im = ax.imshow(heatmap.values.astype(float), aspect="auto", cmap="RdBu_r")
    ax.set_xticks(range(len(stage_cols)))
    ax.set_xticklabels(stage_cols, rotation=90, fontsize=7)
    ax.set_yticks(range(len(cell_types)))
    ax.set_yticklabels(cell_types, fontsize=8)
    plt.colorbar(im, ax=ax, shrink=0.7, label="Mean consensus score")
    ax.set_title("Cell type × Disease stage activation")
    plt.tight_layout()
    out = os.path.join(fig_dir, "celltype_stage_heatmap.pdf")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved: %s", out)


def plot_celltype_program_heatmap(activity_matrix, fig_dir):
    """Heatmap: cell-type disease program cross-scoring."""
    if activity_matrix is None:
        return

    log.info("Plotting cell-type disease program heatmap...")
    fig, ax = plt.subplots(figsize=(max(10, activity_matrix.shape[1] * 0.6),
                                     max(6, activity_matrix.shape[0] * 0.4)))
    vals = activity_matrix.values.astype(float)
    im = ax.imshow(vals, aspect="auto", cmap="RdBu_r")
    ax.set_xticks(range(activity_matrix.shape[1]))
    ax.set_xticklabels([c.replace("ct_", "").replace("_MASLD", "") for c in activity_matrix.columns],
                       rotation=90, fontsize=7)
    ax.set_yticks(range(activity_matrix.shape[0]))
    ax.set_yticklabels(activity_matrix.index, fontsize=8)
    plt.colorbar(im, ax=ax, shrink=0.7, label="Mean consensus score")
    ax.set_title("Cell type × Disease program (cross-scoring)")
    plt.tight_layout()
    out = os.path.join(fig_dir, "celltype_disease_programs.pdf")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved: %s", out)


def plot_pseudotime(adata_sub, dpt_results, fig_dir, celltype_name="Hepatocytes"):
    """UMAP colored by DPT + gene cascade heatmap."""
    if adata_sub is None:
        return

    log.info("Plotting pseudotime for %s...", celltype_name)

    # Recompute UMAP for this subset
    try:
        sc.tl.umap(adata_sub)
    except Exception:
        pass

    if "X_umap" not in adata_sub.obsm:
        log.warning("No UMAP for pseudotime plot")
        return

    dpt = adata_sub.obs["dpt_pseudotime"].values
    valid = ~np.isnan(dpt) & ~np.isinf(dpt)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel 1: UMAP colored by DPT
    ax = axes[0]
    umap = adata_sub.obsm["X_umap"]
    sc_plot = ax.scatter(umap[valid, 0], umap[valid, 1],
                         c=dpt[valid], cmap="magma", s=1, alpha=0.5, rasterized=True)
    plt.colorbar(sc_plot, ax=ax, shrink=0.7, label="DPT")
    ax.set_title(f"{celltype_name} disease pseudotime")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_xticks([])
    ax.set_yticks([])

    # Panel 2: Top pseudotime-correlated genes (heatmap cascade)
    ax = axes[1]
    gene_corr_df = dpt_results.get("gene_correlations")
    if gene_corr_df is not None and len(gene_corr_df) > 0:
        top_genes = gene_corr_df.head(30)["gene"].tolist()
        valid_genes = [g for g in top_genes if g in adata_sub.var_names][:20]

        if valid_genes:
            # Sort cells by DPT
            order = np.argsort(dpt)
            order = order[np.isfinite(dpt[order])]

            X = adata_sub[:, valid_genes].X
            if issparse(X):
                X = X.toarray()
            X_ordered = X[order]

            # Smooth with rolling window for visualization
            window = max(1, len(order) // 100)
            X_smooth = pd.DataFrame(X_ordered).rolling(window, center=True, min_periods=1).mean().values

            # Z-score per gene for heatmap
            X_z = (X_smooth - X_smooth.mean(axis=0)) / (X_smooth.std(axis=0) + 1e-8)

            im = ax.imshow(X_z.T, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
            ax.set_yticks(range(len(valid_genes)))
            ax.set_yticklabels(valid_genes, fontsize=7)
            ax.set_xlabel("Cells ordered by pseudotime →")
            ax.set_title("Gene expression cascade")
            plt.colorbar(im, ax=ax, shrink=0.7, label="Z-score")
        else:
            ax.text(0.5, 0.5, "No genes mapped", transform=ax.transAxes, ha="center")
    else:
        ax.text(0.5, 0.5, "No gene correlations", transform=ax.transAxes, ha="center")

    plt.tight_layout()
    out = os.path.join(fig_dir, f"pseudotime_{celltype_name}.pdf")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved: %s", out)


def plot_cross_method_scatter(aucell_df, scanpy_df, zscore_df, mapped_sigs, fig_dir):
    """Scatter plots comparing AUCell vs scanpy vs z-score for selected signatures."""
    log.info("Plotting cross-method agreement...")

    # Pick first 6 stage signatures for illustration
    stage_sigs = [s for s in mapped_sigs if "vs" in s and not s.startswith("ct_")][:6]
    if not stage_sigs:
        return

    n = len(stage_sigs)
    fig, axes = plt.subplots(2, n, figsize=(4 * n, 8))
    if n == 1:
        axes = axes.reshape(2, 1)

    for i, sig in enumerate(stage_sigs):
        # Row 1: AUCell vs scanpy net
        ax = axes[0, i]
        if sig in aucell_df.columns:
            up_col = f"{sig}_UP"
            dn_col = f"{sig}_DOWN"
            up = scanpy_df[up_col] if up_col in scanpy_df.columns else 0
            dn = scanpy_df[dn_col] if dn_col in scanpy_df.columns else 0
            net_sp = up - dn
            if isinstance(net_sp, pd.Series):
                # Subsample for plotting
                n_plot = min(5000, len(net_sp))
                idx = np.random.choice(len(net_sp), n_plot, replace=False)
                ax.scatter(aucell_df[sig].iloc[idx], net_sp.iloc[idx],
                           s=1, alpha=0.2, rasterized=True)
                ax.set_xlabel("AUCell")
                ax.set_ylabel("score_genes")
        ax.set_title(sig, fontsize=8)

        # Row 2: AUCell vs z-score
        ax = axes[1, i]
        if sig in aucell_df.columns and sig in zscore_df.columns:
            n_plot = min(5000, aucell_df.shape[0])
            idx = np.random.choice(aucell_df.shape[0], n_plot, replace=False)
            ax.scatter(aucell_df[sig].iloc[idx], zscore_df[sig].iloc[idx],
                       s=1, alpha=0.2, rasterized=True)
            ax.set_xlabel("AUCell")
            ax.set_ylabel("z-score")

    plt.tight_layout()
    out = os.path.join(fig_dir, "cross_method_scatter.pdf")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved: %s", out)


def plot_calibration(validation, fig_dir):
    """Pseudobulk score vs condition ordinal."""
    pb = validation.get("pseudobulk_calibration")
    if pb is None:
        return

    log.info("Plotting pseudobulk calibration...")
    # Map to fine condition labels for calibration
    cond_fine_map = {"NAFLD": "MASL", "NASH": "MASH"}
    pb = pb.copy()
    pb["condition"] = pb["condition"].replace(cond_fine_map)
    condition_order = {"Healthy": 0, "MASL": 1, "MASLD": 2, "MASH": 3, "Cirrhotic": 4}
    pb = pb[pb["condition"].isin(condition_order)].copy()
    pb["ordinal"] = pb["condition"].map(condition_order)

    sev_cols = [c for c in pb.columns if "severity" in c]
    if not sev_cols:
        return

    fig, axes = plt.subplots(1, len(sev_cols), figsize=(5 * len(sev_cols), 4))
    if len(sev_cols) == 1:
        axes = [axes]

    for ax, col in zip(axes, sev_cols):
        ax.scatter(pb["ordinal"], pb[col], alpha=0.5, s=20, rasterized=True)
        # Add boxplot overlay
        for ord_val in sorted(pb["ordinal"].unique()):
            vals = pb[pb["ordinal"] == ord_val][col].dropna()
            if len(vals) > 0:
                bp = ax.boxplot([vals], positions=[ord_val], widths=0.3,
                                showfliers=False, patch_artist=True)
                for patch in bp["boxes"]:
                    patch.set_facecolor("lightblue")
                    patch.set_alpha(0.5)

        rho, p = stats.spearmanr(pb["ordinal"], pb[col].fillna(0))
        ax.set_title(f"{col}\nrho={rho:.3f}, p={p:.2e}", fontsize=9)
        ax.set_xticks(list(condition_order.values()))
        ax.set_xticklabels(list(condition_order.keys()), rotation=45, ha="right")
        ax.set_ylabel("Mean severity score")

    plt.tight_layout()
    out = os.path.join(fig_dir, "pseudobulk_calibration.pdf")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved: %s", out)


# ============================================================================
# Main Pipeline
# ============================================================================
def main():
    log.info("=" * 70)
    log.info("Disease Signature Scoring Pipeline")
    log.info("=" * 70)

    # ------------------------------------------------------------------
    # Step 1: Extract signatures from bulk DEGs
    # ------------------------------------------------------------------
    log.info("\n=== STEP 1: Signature extraction ===")

    # 1a. NAS consecutive signatures (k vs k-1)
    nas_consec_path = os.path.join(DISEASE_SIG_DIR, "nas_consecutive_dream.csv")
    log.info("Loading NAS consecutive DEGs: %s", nas_consec_path)
    nas_consec = pd.read_csv(nas_consec_path)
    log.info("  %d rows, contrasts: %s", len(nas_consec), nas_consec["contrast"].unique().tolist())
    nas_consec_sigs = extract_stage_signatures(nas_consec, "contrast")

    # 1b. Fibrosis consecutive signatures (k vs k-1)
    fib_consec_path = os.path.join(DISEASE_SIG_DIR, "fibrosis_consecutive_dream.csv")
    log.info("Loading fibrosis consecutive DEGs: %s", fib_consec_path)
    fib_consec = pd.read_csv(fib_consec_path)
    log.info("  %d rows, contrasts: %s", len(fib_consec), fib_consec["contrast"].unique().tolist())
    fib_consec_sigs = extract_stage_signatures(fib_consec, "contrast")

    # 1c. Cumulative signatures (k vs baseline) for Approach B
    nas_cum_path = os.path.join(DISEASE_SIG_DIR, "nas_score_dream.csv")
    log.info("Loading NAS cumulative DEGs: %s", nas_cum_path)
    nas_cum = pd.read_csv(nas_cum_path)
    nas_cum_sigs = extract_stage_signatures(nas_cum, "contrast")

    fib_cum_path = os.path.join(DISEASE_SIG_DIR, "fibrosis_stage_dream.csv")
    log.info("Loading fibrosis cumulative DEGs: %s", fib_cum_path)
    fib_cum = pd.read_csv(fib_cum_path)
    fib_cum_sigs = extract_stage_signatures(fib_cum, "contrast")

    # 1d. Cell-type-specific signatures from pseudobulk DE
    log.info("Loading cell-type pseudobulk DE signatures...")
    ct_sigs = extract_celltype_signatures(PSEUDOBULK_DIR)

    # Merge all signatures
    all_sigs = {}
    all_sigs.update(nas_consec_sigs)
    all_sigs.update(fib_consec_sigs)
    all_sigs.update(nas_cum_sigs)
    all_sigs.update(fib_cum_sigs)
    all_sigs.update(ct_sigs)
    log.info("Total signature sets: %d", len(all_sigs))

    # Free large DataFrames
    del nas_consec, fib_consec, nas_cum, fib_cum
    gc.collect()

    # ------------------------------------------------------------------
    # Step 2: Load atlas & map gene IDs
    # ------------------------------------------------------------------
    log.info("\n=== STEP 2: Load atlas & map gene IDs ===")
    adata = load_atlas()

    id_to_var = build_gene_id_map(adata)
    mapped_sigs, mapping_stats = map_signatures_to_atlas(all_sigs, id_to_var)

    # Save gene ID mapping and signature stats
    mapping_stats.to_csv(os.path.join(OUT_DIR, "signature_mapping_stats.csv"), index=False)

    # Save signature gene sets
    sig_rows = []
    for sig_name, directions in mapped_sigs.items():
        for direction, genes in directions.items():
            for gene in genes:
                sig_rows.append({"signature": sig_name, "direction": direction, "gene": gene})
    sig_df = pd.DataFrame(sig_rows)
    sig_df.to_csv(os.path.join(OUT_DIR, "signature_gene_sets.csv"), index=False)
    log.info("Saved signature gene sets: %d rows", len(sig_df))

    # Save gene ID mapping
    id_map_rows = [{"versioned_ensembl": f"{k}.x", "unversioned": k, "atlas_var_name": v}
                   for k, v in id_to_var.items()]
    pd.DataFrame(id_map_rows).to_csv(os.path.join(OUT_DIR, "gene_id_mapping.csv"), index=False)

    # Build decoupler nets
    net_combined, net_directional = signatures_to_net(mapped_sigs)
    log.info("Net (combined): %d edges, %d sources", len(net_combined), net_combined["source"].nunique())
    log.info("Net (directional): %d edges, %d sources", len(net_directional), net_directional["source"].nunique())

    # ------------------------------------------------------------------
    # Step 3: Multi-method ensemble scoring
    # ------------------------------------------------------------------
    log.info("\n=== STEP 3: Multi-method ensemble scoring ===")

    # Score methods one at a time, saving + freeing after each to manage RAM

    # AUCell (directional: UP and DOWN scored separately, then netted)
    log.info("--- AUCell ---")
    aucell_df = score_aucell(adata, net_directional, mapped_sigs)
    aucell_df.to_csv(os.path.join(OUT_DIR, "scores_aucell.csv.gz"), compression="gzip")
    log.info("  Saved scores_aucell.csv.gz (%d cols)", aucell_df.shape[1])
    del aucell_df; gc.collect()

    # scanpy score_genes (separate UP and DOWN per signature)
    log.info("--- scanpy score_genes ---")
    scanpy_df = score_scanpy(adata, mapped_sigs)
    scanpy_df.to_csv(os.path.join(OUT_DIR, "scores_scanpy.csv.gz"), compression="gzip")
    log.info("  Saved scores_scanpy.csv.gz (%d cols)", scanpy_df.shape[1])
    del scanpy_df; gc.collect()

    # z-score (combined net with signed weights)
    log.info("--- z-score ---")
    zscore_df = score_zscore(adata, net_combined)
    zscore_df.to_csv(os.path.join(OUT_DIR, "scores_zscore.csv.gz"), compression="gzip")
    log.info("  Saved scores_zscore.csv.gz (%d cols)", zscore_df.shape[1])
    del zscore_df; gc.collect()

    # ------------------------------------------------------------------
    # Step 3b: Pseudobulk-level scoring — all cells (before freeing nets)
    # ------------------------------------------------------------------
    log.info("\n=== STEP 3b: Pseudobulk-level scoring (all cells) ===")
    pb_adata = create_pseudobulk_adata(adata)
    pb_severity = None
    if pb_adata is not None:
        pb_severity = score_pseudobulk(pb_adata, net_combined, net_directional, mapped_sigs)
        del pb_adata; gc.collect()

    # ------------------------------------------------------------------
    # Step 3c: Hepatocyte-only pseudobulk (composition-controlled)
    # ------------------------------------------------------------------
    log.info("\n=== STEP 3c: Hepatocyte-only pseudobulk scoring ===")
    hep_pb = create_pseudobulk_adata(adata, cell_type_filter="Hepatocytes", min_cells=30)
    hep_severity = None
    if hep_pb is not None:
        hep_severity = score_pseudobulk(
            hep_pb, net_combined, net_directional, mapped_sigs,
            prefix="hep_pseudobulk",
        )
        del hep_pb; gc.collect()

    # Free nets
    del net_combined, net_directional; gc.collect()

    # Reload all for consensus (each is ~200-400MB compressed, ~2-5GB in memory)
    log.info("--- Ensemble consensus ---")
    log.info("Reloading per-method scores for consensus...")
    aucell_df = pd.read_csv(os.path.join(OUT_DIR, "scores_aucell.csv.gz"), index_col=0)
    scanpy_df = pd.read_csv(os.path.join(OUT_DIR, "scores_scanpy.csv.gz"), index_col=0)
    zscore_df = pd.read_csv(os.path.join(OUT_DIR, "scores_zscore.csv.gz"), index_col=0)
    consensus_df = compute_ensemble_scores(aucell_df, scanpy_df, zscore_df, mapped_sigs)
    del aucell_df, scanpy_df, zscore_df; gc.collect()

    consensus_df.to_csv(os.path.join(OUT_DIR, "scores_consensus.csv.gz"), compression="gzip")
    log.info("  Saved scores_consensus.csv.gz")

    # ------------------------------------------------------------------
    # Step 4: Composite severity indices
    # ------------------------------------------------------------------
    log.info("\n=== STEP 4: Composite severity indices ===")
    severity_df = compute_severity_indices(consensus_df, mapped_sigs)
    severity_df.to_csv(os.path.join(OUT_DIR, "severity_scores.csv"))
    log.info("Severity scores shape: %s", severity_df.shape)

    # Add severity scores to adata.obs for downstream use
    for col in severity_df.columns:
        adata.obs[col] = severity_df.reindex(adata.obs_names)[col].values

    # ------------------------------------------------------------------
    # Step 5: Cell-type-specific disease programs
    # ------------------------------------------------------------------
    log.info("\n=== STEP 5: Cell-type disease programs ===")
    activity_matrix = compute_celltype_disease_programs(consensus_df, adata)
    if activity_matrix is not None:
        activity_matrix.to_csv(os.path.join(OUT_DIR, "celltype_disease_activity.csv"))

    # ------------------------------------------------------------------
    # Step 6: Disease pseudotime
    # ------------------------------------------------------------------
    log.info("\n=== STEP 6: Disease pseudotime ===")
    pseudotime_results = {}
    for ct in ["Hepatocytes", "Fibroblasts", "Macrophages"]:
        log.info("--- Pseudotime: %s ---", ct)
        adata_sub, dpt_results = compute_pseudotime(adata, severity_df, target_celltype=ct)
        if adata_sub is not None:
            pseudotime_results[ct] = (adata_sub, dpt_results)
            # Save per-cell-type DPT
            dpt_out = adata_sub.obs[["dpt_pseudotime"]].copy()
            for col in severity_df.columns:
                if col in adata_sub.obs.columns:
                    dpt_out[col] = adata_sub.obs[col]
            dpt_out.to_csv(os.path.join(OUT_DIR, f"pseudotime_{ct}.csv"))
            # Save gene correlations
            if dpt_results and dpt_results.get("gene_correlations") is not None:
                dpt_results["gene_correlations"].to_csv(
                    os.path.join(OUT_DIR, f"pseudotime_gene_corr_{ct}.csv"), index=False
                )
            if dpt_results and dpt_results.get("dpt_severity_corr") is not None:
                dpt_results["dpt_severity_corr"].to_csv(
                    os.path.join(OUT_DIR, f"pseudotime_severity_corr_{ct}.csv"), index=False
                )

    # ------------------------------------------------------------------
    # Step 7: Validation
    # ------------------------------------------------------------------
    log.info("\n=== STEP 7: Validation ===")
    # Reload per-method scores for cross-method validation (freed earlier to save RAM)
    log.info("Reloading per-method scores for validation...")
    aucell_df = pd.read_csv(os.path.join(OUT_DIR, "scores_aucell.csv.gz"), index_col=0)
    scanpy_df = pd.read_csv(os.path.join(OUT_DIR, "scores_scanpy.csv.gz"), index_col=0)
    zscore_df = pd.read_csv(os.path.join(OUT_DIR, "scores_zscore.csv.gz"), index_col=0)
    validation = validate_scores(adata, severity_df, consensus_df,
                                 aucell_df, scanpy_df, zscore_df, mapped_sigs)
    del aucell_df, scanpy_df, zscore_df; gc.collect()

    # Save validation outputs
    for key, val in validation.items():
        if isinstance(val, pd.DataFrame):
            val.to_csv(os.path.join(OUT_DIR, f"validation_{key}.csv"), index=True)
            log.info("Saved validation_%s.csv", key)

    # ------------------------------------------------------------------
    # Step 8: Visualization
    # ------------------------------------------------------------------
    log.info("\n=== STEP 8: Visualization ===")
    plot_umap_severity(adata, severity_df, FIG_DIR)
    plot_celltype_violins(adata, severity_df, FIG_DIR)
    plot_celltype_stage_heatmap(consensus_df, adata, FIG_DIR)
    plot_celltype_program_heatmap(activity_matrix, FIG_DIR)
    # Cross-method scatter: skip if per-method scores already freed
    # (validation already computed agreement stats above)
    plot_calibration(validation, FIG_DIR)
    # Pseudobulk-level calibration (expression-based, not cell-level aggregation)
    pb_sev_path = os.path.join(OUT_DIR, "pseudobulk_severity_scores.csv")
    if os.path.exists(pb_sev_path):
        plot_pseudobulk_calibration(pb_sev_path, FIG_DIR)
    # Hepatocyte-only pseudobulk calibration (composition-controlled)
    hep_sev_path = os.path.join(OUT_DIR, "hep_pseudobulk_severity_scores.csv")
    if os.path.exists(hep_sev_path):
        plot_pseudobulk_calibration(hep_sev_path, FIG_DIR)

    for ct, (adata_sub, dpt_results) in pseudotime_results.items():
        plot_pseudotime(adata_sub, dpt_results, FIG_DIR, celltype_name=ct)

    # ------------------------------------------------------------------
    # Save scored atlas (severity scores in .obs)
    # ------------------------------------------------------------------
    log.info("\n=== Saving scored atlas ===")
    scored_path = os.path.join(OUT_DIR, "scored_atlas.h5ad")
    # Only save obs columns (not the full expression matrix) to save space
    adata_out = ad.AnnData(
        obs=adata.obs,
        obsm={"X_umap": adata.obsm["X_umap"]} if "X_umap" in adata.obsm else {},
    )
    adata_out.write_h5ad(scored_path)
    log.info("Saved scored atlas metadata: %s", scored_path)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    log.info("\n" + "=" * 70)
    log.info("PIPELINE COMPLETE")
    log.info("=" * 70)
    log.info("Output directory: %s", OUT_DIR)
    log.info("Signatures extracted: %d total", len(all_sigs))
    log.info("Cells scored: %d", adata.n_obs)
    log.info("Severity columns: %s", list(severity_df.columns))
    if activity_matrix is not None:
        log.info("Cell-type activity matrix: %s", activity_matrix.shape)
    log.info("Pseudotime computed for: %s", list(pseudotime_results.keys()))

    # Verification checklist
    log.info("\n--- Verification ---")
    sig_sizes = mapping_stats.groupby("signature")["n_mapped"].mean()
    low_sigs = sig_sizes[sig_sizes < 70]
    if len(low_sigs) > 0:
        log.warning("Gene sets with <70 mapped genes: %d", len(low_sigs))
    else:
        log.info("All gene sets have >= 70 mapped genes")

    nan_counts = severity_df.isna().sum()
    log.info("NaN counts in severity scores: %s", nan_counts.to_dict())

    concordance = validation.get("concordance")
    if concordance is not None and len(concordance) > 0:
        sig_trends = concordance[concordance["jt_pvalue"] < 0.05]
        log.info("Scores with significant monotonic trend: %d / %d",
                 len(sig_trends), len(concordance))

    agreement = validation.get("cross_method_agreement")
    if agreement is not None and len(agreement) > 0:
        log.info("Cross-method agreement: mean rho=%.3f ± %.3f",
                 agreement["spearman_rho"].mean(), agreement["spearman_rho"].std())

    log.info("Done!")


if __name__ == "__main__":
    main()
