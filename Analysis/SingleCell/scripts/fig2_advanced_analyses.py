#!/usr/bin/env python3
"""
Figure 3 advanced analyses: marker extraction, TF activity (decoupleR), LIANA CCC.
Outputs CSVs to results_gpu_v2/fig2_data/ for downstream R plotting.

Usage:
    sbatch run_fig3_advanced.sh
"""

import os
import sys
import warnings
import logging

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
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
OUT_DIR = os.path.join(RESULTS, "fig2_data")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# GPU init (optional — speeds up h5ad loading via rapids)
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.join(SC_DIR, "scripts"))
try:
    from gpu_utils import init_gpu
    init_gpu()
    log.info("GPU initialized")
except Exception:
    log.info("No GPU or gpu_utils unavailable — running CPU-only")

# ---------------------------------------------------------------------------
# Imports (after GPU init to avoid CUDA ordering issues)
# ---------------------------------------------------------------------------
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

# ============================================================================
# 1. Load h5ad (with fix for missing 'ordered' attribute on categoricals)
# ============================================================================
log.info("Loading h5ad: %s", H5AD)

# Fix missing 'ordered' attribute before reading with scanpy
import h5py

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

fix_ordered_attrs(H5AD)

adata = sc.read_h5ad(H5AD)
log.info("Loaded %d cells x %d genes", adata.n_obs, adata.n_vars)

# Ensure raw counts or normalized data available
# Use .raw if present (often stores log-normalized data)
if adata.raw is not None:
    log.info("Using adata.raw for expression values")
    adata_expr = adata.raw.to_adata()
else:
    adata_expr = adata

# ============================================================================
# 2. Marker gene dotplot data extraction
# ============================================================================
log.info("=== Step 1: Extracting marker gene dotplot data ===")

MARKERS = {
    "Hepatocytes": ["ALB", "APOA1"],
    "Cholangiocytes": ["KRT19", "SOX9"],
    "Endothelial cells": ["PECAM1", "CDH5"],
    "Fibroblasts": ["COL1A1", "ACTA2"],
    "Macrophages": ["CD68", "MARCO"],
    "T cells": ["CD3D"],
    "B cells": ["CD79A"],
    "Resident NK": ["GNLY"],
    "Plasma cells": ["JCHAIN"],
}

marker_genes = []
for genes in MARKERS.values():
    marker_genes.extend(genes)
marker_genes = list(dict.fromkeys(marker_genes))  # deduplicate, preserve order

# Find which marker genes are in the dataset (may use var_names or gene_ids)
var_names = adata_expr.var_names.tolist()

# Try direct match first
found_markers = [g for g in marker_genes if g in var_names]

# Check for a symbol column in .var
symbol_col = None
for col in ["gene_name", "gene_symbols", "gene_short_name", "feature_name"]:
    if col in adata_expr.var.columns:
        symbol_col = col
        break

# If few found by direct match, try matching via symbol column
if len(found_markers) < len(marker_genes) // 2 and symbol_col:
    log.info("Using var column '%s' for gene symbol matching", symbol_col)
    symbol_map = dict(zip(adata_expr.var[symbol_col], adata_expr.var_names))
    found_markers = []
    for g in marker_genes:
        if g in var_names:
            found_markers.append(g)
        elif g in symbol_map:
            found_markers.append(symbol_map[g])

log.info("Found %d / %d marker genes in dataset", len(found_markers), len(marker_genes))

# Build reverse map: var_name -> display symbol
if symbol_col and symbol_col in adata_expr.var.columns:
    var_to_symbol = dict(zip(adata_expr.var_names, adata_expr.var[symbol_col]))
else:
    var_to_symbol = {g: g for g in var_names}

# Compute mean expression and pct expressing per cell type
cell_type_col = "cell_type"
cell_types = adata_expr.obs[cell_type_col].unique().tolist()

dotplot_rows = []
for ct in cell_types:
    mask = adata_expr.obs[cell_type_col] == ct
    n_cells = mask.sum()
    if n_cells == 0:
        continue
    subset = adata_expr[mask, :]
    for gene_var in found_markers:
        if gene_var not in adata_expr.var_names:
            continue
        gene_idx = adata_expr.var_names.get_loc(gene_var)
        expr = subset[:, gene_idx].X
        if hasattr(expr, "toarray"):
            expr = expr.toarray().flatten()
        else:
            expr = np.asarray(expr).flatten()
        mean_expr = float(np.mean(expr))
        pct_expr = float(np.sum(expr > 0) / n_cells * 100)
        symbol = var_to_symbol.get(gene_var, gene_var)
        dotplot_rows.append({
            "cell_type": ct,
            "gene": symbol,
            "gene_var": gene_var,
            "mean_expr": mean_expr,
            "pct_expressing": pct_expr,
            "n_cells": int(n_cells),
        })

dotplot_df = pd.DataFrame(dotplot_rows)
out_marker = os.path.join(OUT_DIR, "marker_dotplot_data.csv")
dotplot_df.to_csv(out_marker, index=False)
log.info("Marker dotplot data saved: %s (%d rows)", out_marker, len(dotplot_df))

# ============================================================================
# 3. TF regulon activity via decoupleR
# ============================================================================
log.info("=== Step 2: TF activity via decoupleR ===")

import decoupler as dc

# Get CollecTRI network (curated TF-target interactions) — decoupler 2.x API
log.info("Fetching CollecTRI TF-target network...")
net = dc.op.collectri(organism="human", remove_complexes=False)
log.info("CollecTRI network: %d interactions, %d TFs", len(net), net["source"].nunique())

# Run multivariate linear model (MLM) on full cells
# MLM is fast: matrix multiplication, scales linearly with cells
# Subsample to ~200K cells for TF activity (stratified by cell_type + condition)
# to avoid OOM on 1.17M cells
MAX_CELLS_TF = 200000
if adata_expr.n_obs > MAX_CELLS_TF:
    log.info("Subsampling from %d to %d cells for TF activity...", adata_expr.n_obs, MAX_CELLS_TF)
    sc.pp.subsample(adata_expr, n_obs=MAX_CELLS_TF, random_state=42, copy=False)
    log.info("Subsampled to %d cells", adata_expr.n_obs)

log.info("Running decoupleR MLM on %d cells...", adata_expr.n_obs)
dc.mt.mlm(
    data=adata_expr,
    net=net,
    verbose=True,
)
log.info("MLM complete. Result stored in adata_expr.obsm")

# Extract TF activity scores (decoupler 2.x uses score_mlm / padj_mlm keys)
tf_key = "score_mlm" if "score_mlm" in adata_expr.obsm else "mlm_estimate"
pv_key = "padj_mlm" if "padj_mlm" in adata_expr.obsm else "mlm_pvals"
log.info("TF obsm keys: %s", list(adata_expr.obsm.keys()))
tf_acts = adata_expr.obsm[tf_key]
tf_pvals = adata_expr.obsm[pv_key]

# Build condition column — use adata_expr (which may be subsampled)
# Transfer metadata from adata to adata_expr for the subsampled cells
adata_expr.obs["cell_type"] = adata.obs.loc[adata_expr.obs_names, "cell_type"].values
if "condition_harmonized" in adata.obs.columns:
    conditions = adata.obs.loc[adata_expr.obs_names, "condition_harmonized"].values
elif "condition" in adata.obs.columns:
    conditions = adata.obs.loc[adata_expr.obs_names, "condition"].values
else:
    conditions = np.array(["Unknown"] * adata_expr.n_obs)

# Map to MASLD/Control
condition_map = {"Healthy": "Control", "MASLD": "MASLD", "NAFLD": "MASLD", "NASH": "MASLD", "Mixed": "Control"}
disease_status = pd.Series([condition_map.get(str(c), "Unknown") for c in conditions],
                           index=adata_expr.obs.index)

# Aggregate TF activity per cell_type x condition
log.info("Aggregating TF activity per cell type x condition...")
cell_types_obs = adata_expr.obs[cell_type_col].values
tf_names = tf_acts.columns.tolist() if hasattr(tf_acts, "columns") else [f"TF_{i}" for i in range(tf_acts.shape[1])]

# Convert to DataFrame if needed
if not isinstance(tf_acts, pd.DataFrame):
    tf_acts_df = pd.DataFrame(tf_acts, index=adata_expr.obs.index, columns=tf_names)
else:
    tf_acts_df = tf_acts

tf_acts_df["cell_type"] = cell_types_obs
tf_acts_df["condition"] = disease_status.values

# Filter to Control and MASLD only
tf_acts_df = tf_acts_df[tf_acts_df["condition"].isin(["Control", "MASLD"])]

# Compute mean activity and do t-test per cell_type x TF
from scipy.stats import ttest_ind

tf_results = []
for ct in tf_acts_df["cell_type"].unique():
    ct_data = tf_acts_df[tf_acts_df["cell_type"] == ct]
    ctrl = ct_data[ct_data["condition"] == "Control"]
    masld = ct_data[ct_data["condition"] == "MASLD"]
    if len(ctrl) < 3 or len(masld) < 3:
        continue
    for tf in tf_names:
        ctrl_vals = ctrl[tf].values
        masld_vals = masld[tf].values
        mean_ctrl = float(np.nanmean(ctrl_vals))
        mean_masld = float(np.nanmean(masld_vals))
        try:
            stat, pval = ttest_ind(masld_vals, ctrl_vals, equal_var=False, nan_policy="omit")
        except Exception:
            stat, pval = np.nan, np.nan
        tf_results.append({
            "TF": tf,
            "cell_type": ct,
            "mean_activity_control": mean_ctrl,
            "mean_activity_masld": mean_masld,
            "activity_diff": mean_masld - mean_ctrl,
            "tstat": float(stat) if not np.isnan(stat) else np.nan,
            "pvalue": float(pval) if not np.isnan(pval) else np.nan,
        })

tf_df = pd.DataFrame(tf_results)

# FDR correction per cell type
from statsmodels.stats.multitest import multipletests

for ct in tf_df["cell_type"].unique():
    mask = tf_df["cell_type"] == ct
    pvals = tf_df.loc[mask, "pvalue"].values
    valid = ~np.isnan(pvals)
    if valid.sum() > 0:
        _, padj, _, _ = multipletests(pvals[valid], method="fdr_bh")
        padj_full = np.full(len(pvals), np.nan)
        padj_full[valid] = padj
        tf_df.loc[mask, "padj"] = padj_full

out_tf = os.path.join(OUT_DIR, "tf_activity_per_celltype_condition.csv")
tf_df.to_csv(out_tf, index=False)
log.info("TF activity saved: %s (%d rows, %d TFs)", out_tf, len(tf_df), tf_df["TF"].nunique())

# ============================================================================
# 4. LIANA cell-cell communication
# ============================================================================
log.info("=== Step 3: LIANA cell-cell communication ===")

import liana as li

# LIANA needs log-normalized data in .X and raw counts somewhere
# Subsample to manageable size for LIANA (200K cells)
MAX_CELLS_LIANA = 200000
if adata.n_obs > MAX_CELLS_LIANA:
    log.info("Subsampling from %d to %d cells for LIANA...", adata.n_obs, MAX_CELLS_LIANA)
    sc.pp.subsample(adata, n_obs=MAX_CELLS_LIANA, random_state=42, copy=False)
adata_liana = adata.copy()

# Add disease_status to obs
adata_liana.obs["disease_status"] = disease_status.values

# Filter to Control and MASLD
adata_liana = adata_liana[adata_liana.obs["disease_status"].isin(["Control", "MASLD"])].copy()

# LIANA needs cell_type as a categorical
adata_liana.obs[cell_type_col] = adata_liana.obs[cell_type_col].astype("category")

# Run LIANA separately for MASLD and Control
log.info("Running LIANA for MASLD cells (%d)...",
         (adata_liana.obs["disease_status"] == "MASLD").sum())

liana_results = {}
for cond in ["MASLD", "Control"]:
    adata_cond = adata_liana[adata_liana.obs["disease_status"] == cond].copy()
    n_cells = adata_cond.n_obs
    log.info("LIANA %s: %d cells", cond, n_cells)

    # Filter cell types with too few cells
    ct_counts = adata_cond.obs[cell_type_col].value_counts()
    valid_cts = ct_counts[ct_counts >= 10].index.tolist()
    adata_cond = adata_cond[adata_cond.obs[cell_type_col].isin(valid_cts)].copy()
    adata_cond.obs[cell_type_col] = adata_cond.obs[cell_type_col].cat.remove_unused_categories()

    log.info("LIANA %s: %d cells after filtering to %d cell types",
             cond, adata_cond.n_obs, len(valid_cts))

    try:
        li.mt.rank_aggregate(
            adata_cond,
            groupby=cell_type_col,
            resource_name="consensus",
            expr_prop=0.1,
            verbose=True,
            use_raw=True if adata_cond.raw is not None else False,
        )
        res = adata_cond.uns["liana_res"].copy()
        res["condition"] = cond
        liana_results[cond] = res
        log.info("LIANA %s: %d interactions", cond, len(res))
    except Exception as e:
        log.warning("LIANA %s failed: %s", cond, str(e))
        liana_results[cond] = None

# Merge and compute differential interactions
if liana_results.get("MASLD") is not None and liana_results.get("Control") is not None:
    masld_res = liana_results["MASLD"]
    ctrl_res = liana_results["Control"]

    # Save combined raw results
    combined = pd.concat([masld_res, ctrl_res], ignore_index=True)
    out_liana_raw = os.path.join(OUT_DIR, "liana_raw_interactions.csv")
    combined.to_csv(out_liana_raw, index=False)
    log.info("LIANA raw interactions saved: %s", out_liana_raw)

    # Compute differential: merge on source-target-ligand-receptor
    merge_keys = ["source", "target", "ligand_complex", "receptor_complex"]

    # Use magnitude_rank as the primary score (lower = stronger interaction)
    score_col = "magnitude_rank" if "magnitude_rank" in masld_res.columns else "liana_rank"

    masld_slim = masld_res[merge_keys + [score_col]].rename(columns={score_col: "score_masld"})
    ctrl_slim = ctrl_res[merge_keys + [score_col]].rename(columns={score_col: "score_control"})

    diff = masld_slim.merge(ctrl_slim, on=merge_keys, how="outer")
    diff["score_masld"] = diff["score_masld"].fillna(1.0)
    diff["score_control"] = diff["score_control"].fillna(1.0)
    # CONVENTION (A5 fix, 2026-06-20): magnitude_rank is LOWER = stronger, so an
    # interaction stronger in MASLD has score_masld < score_control. We define
    # score_diff = score_control - score_masld so that POSITIVE = stronger/enriched
    # in MASLD and NEGATIVE = stronger in Control. All consumers (R + Python) must
    # use "score_diff > 0 == MASLD-enriched". (Previously score_masld - score_control,
    # which inverted the intuitive sign and was read backwards by the R consumers.)
    diff["score_diff"] = diff["score_control"] - diff["score_masld"]

    # --- Regression guard (formula-based, pair-agnostic): by construction a pair
    # stronger in MASLD has a LOWER magnitude_rank in MASLD (score_masld <
    # score_control) and MUST therefore yield score_diff > 0 under this convention.
    # This pins the formula sign without relying on any single biological pair:
    # individual canonical ligands are unreliable anchors here because magnitude_rank
    # is normalized WITHIN each condition, so secreted inflammatory ligands
    # (SPP1/TIMP1/IL1B) wash out and only ECM-integrin axes robustly rank MASLD-up.
    # A future re-inversion to score_masld - score_control would flip these and trip
    # this guard.
    _masld_stronger = diff["score_masld"] < diff["score_control"]
    if _masld_stronger.any() and not bool((diff.loc[_masld_stronger, "score_diff"] > 0).all()):
        raise AssertionError(
            "score_diff convention inverted: pairs with lower MASLD magnitude_rank "
            "(score_masld < score_control) must have score_diff > 0 (MASLD-enriched)."
        )
    log.info("A5 sign check OK: all %d MASLD-stronger pairs have score_diff > 0",
             int(_masld_stronger.sum()))

    out_liana_diff = os.path.join(OUT_DIR, "liana_differential_interactions.csv")
    diff.to_csv(out_liana_diff, index=False)
    log.info("LIANA differential interactions saved: %s (%d pairs)", out_liana_diff, len(diff))

    # Aggregate communication strength per source-target pair
    agg_rows = []
    for (src, tgt), grp in diff.groupby(["source", "target"]):
        agg_rows.append({
            "source": src,
            "target": tgt,
            "n_interactions": len(grp),
            "mean_score_masld": float(grp["score_masld"].mean()),
            "mean_score_control": float(grp["score_control"].mean()),
            # Convention: score_diff = score_control - score_masld, so > 0 = MASLD-enriched.
            "n_stronger_masld": int((grp["score_diff"] > 0.1).sum()),
            "n_stronger_control": int((grp["score_diff"] < -0.1).sum()),
        })
    agg_df = pd.DataFrame(agg_rows)
    out_liana_agg = os.path.join(OUT_DIR, "liana_aggregated_communication.csv")
    agg_df.to_csv(out_liana_agg, index=False)
    log.info("LIANA aggregated communication saved: %s", out_liana_agg)

else:
    log.warning("LIANA differential analysis skipped (missing results)")
    # Write empty files so downstream R script doesn't fail
    for fname in ["liana_differential_interactions.csv", "liana_aggregated_communication.csv"]:
        pd.DataFrame().to_csv(os.path.join(OUT_DIR, fname), index=False)

# ============================================================================
# Done
# ============================================================================
log.info("=== All fig3 advanced analyses complete ===")
log.info("Output directory: %s", OUT_DIR)
for f in sorted(os.listdir(OUT_DIR)):
    fpath = os.path.join(OUT_DIR, f)
    size_mb = os.path.getsize(fpath) / 1e6
    log.info("  %s (%.1f MB)", f, size_mb)
