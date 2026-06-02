#!/usr/bin/env python3
"""
elatus_lncrna_profiling.py
Module 3: Cell-Type-Specific ncRNA Profiling (ELATUS-inspired)

Applies an ELATUS-inspired dropout correction to lncRNAs in the existing
scRNA-seq atlas, then profiles cell-type-specific expression, computes
tau specificity indices, and runs pseudobulk DE (Wilcoxon) for lncRNAs.

NOTE: This is an approximation of the full ELATUS approach (Bonder et al.
2024 Nat Genet), which re-quantifies BAMs with intron-aware references.
Here we apply a per-cell-type multiplicative correction factor derived
from the ratio of expected bulk lncRNA expression to observed scRNA-seq
lncRNA expression. This corrects for systematic lncRNA dropout without
requiring BAM reprocessing.

Inputs:
  - Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad
  - GENCODE v49 GTF (gzipped)
  - Dream bulk results (for expected lncRNA fraction)

Outputs (all to Analysis/SingleCell/results_gpu_v2/elatus/):
  - lncrna_celltype_expression.csv
  - lncrna_celltype_specific.csv
  - lncrna_celltype_de.csv
  - elatus_qc_metrics.csv
"""

import argparse
import gzip
import logging
import os
import sys
import time
import warnings
from datetime import datetime

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

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
DEFAULT_ATLAS = os.path.join(RESULTS, "integrated_atlas.h5ad")
DEFAULT_OUT = os.path.join(RESULTS, "elatus")
GENCODE_GTF = (
    "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/"
    "gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
)
DREAM_RESULTS = os.path.join(
    BASE,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/"
    "integration/dream_results.csv",
)
SCENIC_REGULONS = os.path.join(
    BASE,
    "Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv",
)

# Known column name candidates
CELLTYPE_COLS = [
    "celltypist_label", "cell_type", "celltype",
    "cell_type_fine", "leiden", "louvain",
]
CONDITION_COLS = [
    "condition_harmonized", "condition", "disease",
    "diagnosis", "status", "group",
]
SAMPLE_COLS = ["sample", "sample_id", "donor", "patient", "orig.ident"]


# ============================================================================
# Section 1 — Parse GENCODE GTF for ncRNA gene biotypes
# ============================================================================
def parse_gencode_biotypes(gtf_path):
    """Parse GENCODE GTF to extract gene name → biotype mapping.

    Returns a dict of {biotype: set(gene_names)} for ncRNA categories.
    """
    log.info("=" * 70)
    log.info("SECTION 1: Parsing GENCODE GTF for ncRNA biotypes")
    log.info("=" * 70)
    t0 = time.time()

    gene_biotype = {}  # gene_name -> gene_type
    ncrna_types = {"lncRNA", "miRNA", "snoRNA", "snRNA"}

    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = fields[8]
            # Extract gene_name and gene_type from attributes
            gn = _extract_attr(attrs, "gene_name")
            gt = _extract_attr(attrs, "gene_type")
            gid = _extract_attr(attrs, "gene_id")
            if gn and gt:
                gene_biotype[gn] = gt
            # Also store by gene_id (without version) for Ensembl-based atlases
            if gid and gt:
                gene_biotype[gid.split(".")[0]] = gt

    # Build per-biotype sets
    biotype_sets = {}
    for bt in ncrna_types:
        genes = {g for g, t in gene_biotype.items() if t == bt}
        biotype_sets[bt] = genes
        log.info("  %s: %d genes", bt, len(genes))

    elapsed = time.time() - t0
    log.info("GTF parsing complete in %.1f s — %d total genes", elapsed, len(gene_biotype))
    return gene_biotype, biotype_sets


def _extract_attr(attrs_str, key):
    """Extract a single attribute value from a GTF attributes string."""
    # Format: key "value"; key "value"; ...
    needle = key + ' "'
    idx = attrs_str.find(needle)
    if idx == -1:
        return None
    start = idx + len(needle)
    end = attrs_str.index('"', start)
    return attrs_str[start:end]


# ============================================================================
# Section 2 — Load scRNA-seq atlas and identify lncRNA genes
# ============================================================================
def load_atlas(atlas_path):
    """Load h5ad atlas and identify metadata columns."""
    log.info("=" * 70)
    log.info("SECTION 2: Loading scRNA-seq atlas")
    log.info("=" * 70)
    t0 = time.time()

    import scanpy as sc
    import h5py

    # Fix missing 'ordered' attribute on categoricals (known issue)
    with h5py.File(atlas_path, "a") as f:
        for group_name in ["obs", "var"]:
            if group_name not in f:
                continue
            grp = f[group_name]
            for key in grp.keys():
                item = grp[key]
                if isinstance(item, h5py.Group) and "categories" in item:
                    if "ordered" not in item.attrs:
                        item.attrs["ordered"] = False

    adata = sc.read_h5ad(atlas_path)
    log.info("Loaded %d cells x %d genes", adata.n_obs, adata.n_vars)
    log.info("obs columns: %s", list(adata.obs.columns))
    log.info("var columns: %s", list(adata.var.columns))
    if adata.uns:
        log.info("uns keys: %s", list(adata.uns.keys()))

    # Identify cell type column
    ct_col = _find_column(adata.obs.columns, CELLTYPE_COLS)
    log.info("Cell type column: %s", ct_col)
    if ct_col:
        cts = adata.obs[ct_col].value_counts()
        log.info("Cell types (%d):\n%s", len(cts), cts.to_string())

    # Identify condition column
    cond_col = _find_column(adata.obs.columns, CONDITION_COLS)
    log.info("Condition column: %s", cond_col)
    if cond_col:
        log.info("Conditions: %s", adata.obs[cond_col].value_counts().to_dict())

    # Identify sample column
    sample_col = _find_column(adata.obs.columns, SAMPLE_COLS)
    log.info("Sample column: %s", sample_col)

    elapsed = time.time() - t0
    log.info("Atlas loaded in %.1f s", elapsed)
    return adata, ct_col, cond_col, sample_col


def _find_column(columns, candidates):
    """Return the first matching column name from candidates."""
    for c in candidates:
        if c in columns:
            return c
    return None


def identify_lncrna_genes(adata, lncrna_set):
    """Identify lncRNA genes present in the atlas.

    Tries matching on var_names directly, then on common symbol columns.
    Returns a boolean mask over adata.var and the matched gene names.
    """
    var_names = set(adata.var_names)
    direct = var_names & lncrna_set
    log.info("lncRNA genes in atlas (direct var_names match): %d / %d GENCODE lncRNAs",
             len(direct), len(lncrna_set))

    # Also check gene symbol columns in .var
    symbol_col = None
    for col in ["gene_name", "gene_symbols", "gene_short_name", "feature_name"]:
        if col in adata.var.columns:
            symbol_col = col
            break

    if symbol_col and len(direct) < 100:
        symbols = set(adata.var[symbol_col].dropna().astype(str))
        sym_match = symbols & lncrna_set
        log.info("lncRNA genes via .var['%s']: %d", symbol_col, len(sym_match))
        if len(sym_match) > len(direct):
            # Build mask using symbol column
            mask = adata.var[symbol_col].isin(lncrna_set)
            matched_names = set(adata.var.loc[mask, symbol_col])
            log.info("Using symbol column '%s' for lncRNA identification", symbol_col)
            return mask, matched_names, symbol_col

    # Use direct var_names match
    mask = pd.Series(adata.var_names).isin(lncrna_set).values
    return mask, direct, None


# ============================================================================
# Section 3 — ELATUS-inspired dropout correction
# ============================================================================
def elatus_correction(adata, lncrna_mask, ct_col, dream_path):
    """Extract raw lncRNA expression matrix and compute QC metrics.

    NOTE (C2 fix): ELATUS multiplicative correction has been REMOVED.
    The original correction factors were inverted (0.03-0.10 instead of >1.0),
    reducing lncRNA counts 90-97% instead of amplifying them. Root cause:
    bulk AveExpr-derived lncRNA fraction (~0.3%) vs scRNA UMI fractions
    (~4-11%) are fundamentally different quantities. A proper ELATUS correction
    requires BAM reprocessing with intron-aware references (Bonder et al. 2024).
    Using uncorrected raw counts for tau and DE is more conservative and honest.

    Returns raw lncRNA expression matrix (cells x lncRNA genes) and
    a QC metrics DataFrame with CF=1.0 (no correction applied).
    """
    log.info("=" * 70)
    log.info("SECTION 3: lncRNA expression extraction (ELATUS correction DISABLED)")
    log.info("=" * 70)
    t0 = time.time()

    # Get raw counts if available
    X = _get_raw_counts(adata)
    log.info("Expression matrix: %s, dtype=%s", type(X).__name__,
             X.dtype if hasattr(X, "dtype") else "unknown")

    lncrna_idx = np.where(lncrna_mask)[0]

    # Per-cell-type QC metrics (no correction applied)
    cell_types = adata.obs[ct_col].unique()
    qc_records = []

    for ct in cell_types:
        ct_mask = (adata.obs[ct_col] == ct).values
        n_cells = ct_mask.sum()
        if n_cells < 10:
            log.info("  %s: skipping (%d cells)", ct, n_cells)
            continue

        X_ct = X[ct_mask, :]

        # Total UMIs per cell
        if sparse.issparse(X_ct):
            total_umi = np.array(X_ct.sum(axis=1)).ravel()
            lncrna_umi = np.array(X_ct[:, lncrna_idx].sum(axis=1)).ravel()
        else:
            total_umi = X_ct.sum(axis=1)
            lncrna_umi = X_ct[:, lncrna_idx].sum(axis=1)

        total_sum = total_umi.sum()
        lncrna_sum = lncrna_umi.sum()
        observed_frac = lncrna_sum / total_sum if total_sum > 0 else 0.0

        qc_records.append({
            "cell_type": ct,
            "n_cells": n_cells,
            "total_umi": int(total_sum),
            "lncrna_umi": int(lncrna_sum),
            "lncrna_frac_observed": round(observed_frac, 6),
            "correction_factor": 1.0,
            "correction_applied": False,
        })
        log.info("  %s: n=%d, lncRNA_frac=%.4f, CF=1.0 (no correction)",
                 ct, n_cells, observed_frac)

    qc_df = pd.DataFrame(qc_records)

    # Extract raw lncRNA matrix (no correction)
    log.info("Extracting raw lncRNA counts (no ELATUS correction)...")
    X_lncrna = X[:, lncrna_idx]
    if sparse.issparse(X_lncrna):
        X_lncrna = X_lncrna.copy().astype(np.float64)
    else:
        X_lncrna = X_lncrna.astype(np.float64)

    elapsed = time.time() - t0
    log.info("lncRNA extraction complete in %.1f s", elapsed)
    return X_lncrna, qc_df, {}


def _get_raw_counts(adata):
    """Retrieve the best available count matrix (prefer raw/counts layer)."""
    if "counts" in adata.layers:
        log.info("Using adata.layers['counts']")
        return adata.layers["counts"]
    if adata.raw is not None:
        log.info("Using adata.raw.X")
        return adata.raw.X
    log.info("Using adata.X (no raw/counts layer found)")
    return adata.X


def _compute_bulk_lncrna_fraction(dream_path, lncrna_names):
    """Estimate expected lncRNA expression fraction from dream bulk results.

    Uses AveExpr (mean log2-CPM) to estimate the fraction of total expression
    attributable to lncRNAs in bulk RNA-seq.
    """
    if not os.path.exists(dream_path):
        log.warning("Dream results not found at %s — using default fraction 0.05", dream_path)
        return 0.05

    try:
        dream = pd.read_csv(dream_path, usecols=["gene", "AveExpr"])
        # Strip version numbers from gene IDs if present
        dream["gene_base"] = dream["gene"].str.split(".").str[0]

        # Convert AveExpr (log2-CPM) to linear scale for fraction estimation
        dream["linear_expr"] = 2.0 ** dream["AveExpr"]

        is_lncrna = dream["gene_base"].isin(lncrna_names) | dream["gene"].isin(lncrna_names)
        lncrna_expr = dream.loc[is_lncrna, "linear_expr"].sum()
        total_expr = dream["linear_expr"].sum()
        frac = lncrna_expr / total_expr if total_expr > 0 else 0.05

        n_matched = is_lncrna.sum()
        log.info("Dream: %d/%d genes matched as lncRNA, bulk fraction=%.4f",
                 n_matched, len(dream), frac)
        return max(frac, 0.001)  # floor at 0.1% to avoid division issues
    except Exception as e:
        log.warning("Failed to parse dream results: %s — using default 0.05", e)
        return 0.05


# ============================================================================
# Section 4 — Cell-type expression profiling
# ============================================================================
def celltype_expression(adata, X_lncrna, lncrna_genes, ct_col, cond_col):
    """Compute mean expression and detection rate per cell type per condition.

    Returns a long-format DataFrame.
    """
    log.info("=" * 70)
    log.info("SECTION 4: Cell-type expression profiling")
    log.info("=" * 70)
    t0 = time.time()

    records = []
    cell_types = adata.obs[ct_col].unique()

    if cond_col:
        conditions = adata.obs[cond_col].unique()
    else:
        conditions = ["all"]

    for ct in cell_types:
        ct_mask = (adata.obs[ct_col] == ct).values
        if ct_mask.sum() < 10:
            continue

        for cond in conditions:
            if cond_col and cond != "all":
                mask = ct_mask & (adata.obs[cond_col] == cond).values
            else:
                mask = ct_mask

            n_cells = mask.sum()
            if n_cells < 5:
                continue

            expr = X_lncrna[mask, :]
            if sparse.issparse(expr):
                det_rate = np.asarray((expr > 0).mean(axis=0)).ravel()
                mean_expr = np.asarray(expr.mean(axis=0)).ravel()
            else:
                det_rate = (expr > 0).mean(axis=0)
                mean_expr = expr.mean(axis=0)

            for i, gene in enumerate(lncrna_genes):
                records.append({
                    "gene": gene,
                    "cell_type": ct,
                    "condition": cond,
                    "detection_rate": round(float(det_rate[i]), 6),
                    "mean_expr": round(float(mean_expr[i]), 6),
                    "n_cells": n_cells,
                })

    df = pd.DataFrame(records)
    elapsed = time.time() - t0
    log.info("Expression profiling: %d rows in %.1f s", len(df), elapsed)
    return df


# ============================================================================
# Section 5 — Tau specificity index
# ============================================================================
def compute_tau(adata, X_lncrna, lncrna_genes, ct_col):
    """Compute tau tissue-specificity index per lncRNA across cell types.

    tau = sum(1 - x_i/max(x)) / (n - 1)
    where x_i is mean expression in cell type i.

    tau > 0.8: cell-type-specific
    tau < 0.3: ubiquitous
    """
    log.info("=" * 70)
    log.info("SECTION 5: Tau specificity index")
    log.info("=" * 70)
    t0 = time.time()

    cell_types = sorted(adata.obs[ct_col].unique())
    # Filter to cell types with sufficient cells
    cell_types = [ct for ct in cell_types
                  if (adata.obs[ct_col] == ct).sum() >= 10]
    n_ct = len(cell_types)
    log.info("Computing tau across %d cell types", n_ct)

    if n_ct < 2:
        log.warning("Need >= 2 cell types for tau; skipping")
        return pd.DataFrame()

    # Build mean expression matrix (genes x cell_types)
    mean_mat = np.zeros((len(lncrna_genes), n_ct))
    for j, ct in enumerate(cell_types):
        ct_mask = (adata.obs[ct_col] == ct).values
        ct_mean = X_lncrna[ct_mask, :].mean(axis=0)
        if sparse.issparse(X_lncrna):
            mean_mat[:, j] = np.asarray(ct_mean).ravel()
        else:
            mean_mat[:, j] = ct_mean

    records = []
    for i, gene in enumerate(lncrna_genes):
        x = mean_mat[i, :]
        x_max = x.max()
        if x_max == 0:
            tau = 0.0
            best_ct = cell_types[0]
            best_expr = 0.0
            n_expressed = 0
        else:
            hat = x / x_max
            tau = float(np.sum(1.0 - hat) / (n_ct - 1))
            best_idx = int(np.argmax(x))
            best_ct = cell_types[best_idx]
            best_expr = float(x_max)
            n_expressed = int((x > 0).sum())

        records.append({
            "gene": gene,
            "tau": round(tau, 4),
            "most_specific_celltype": best_ct,
            "most_specific_expr": round(best_expr, 6),
            "n_celltypes_expressed": n_expressed,
        })

    df = pd.DataFrame(records)
    n_specific = (df["tau"] > 0.8).sum()
    n_ubiq = (df["tau"] < 0.3).sum()
    log.info("Tau: %d specific (>0.8), %d ubiquitous (<0.3), %d intermediate",
             n_specific, n_ubiq, len(df) - n_specific - n_ubiq)

    elapsed = time.time() - t0
    log.info("Tau computation complete in %.1f s", elapsed)
    return df


# ============================================================================
# Section 6 — Pseudobulk DE for lncRNAs (Wilcoxon)
# ============================================================================
def pseudobulk_de(adata, X_lncrna, lncrna_genes, ct_col, cond_col, sample_col):
    """Pseudobulk differential expression for lncRNAs per cell type.

    Aggregates counts to pseudobulk (sum per sample per cell type), then
    runs Wilcoxon rank-sum (Mann-Whitney U) between disease and control
    pseudobulk samples. Full limma-voom pseudobulk DE requires R; this
    provides a nonparametric approximation.
    """
    log.info("=" * 70)
    log.info("SECTION 6: Pseudobulk DE for lncRNAs")
    log.info("=" * 70)
    t0 = time.time()

    if not cond_col:
        log.warning("No condition column found — skipping pseudobulk DE")
        return pd.DataFrame()
    if not sample_col:
        log.warning("No sample column found — skipping pseudobulk DE")
        return pd.DataFrame()

    # Identify disease vs control conditions
    cond_values = adata.obs[cond_col].value_counts()
    log.info("Condition distribution:\n%s", cond_values.to_string())

    control_labels = {"Healthy", "Control", "Normal", "healthy", "control", "normal"}
    disease_labels = {"MASLD", "NASH", "NAFL", "Cirrhotic", "Disease", "Fibrosis",
                      "masld", "nash", "nafl", "disease", "cirrhotic", "fibrosis"}

    ctrl = [c for c in cond_values.index if c in control_labels]
    dis = [c for c in cond_values.index if c in disease_labels]

    if not ctrl or not dis:
        # Fallback: treat anything not "Healthy"/"Control" as disease
        ctrl = [c for c in cond_values.index
                if any(h in str(c).lower() for h in ["healthy", "control", "normal"])]
        dis = [c for c in cond_values.index if c not in ctrl]

    if not ctrl or not dis:
        log.warning("Cannot determine disease vs control groups — skipping DE")
        return pd.DataFrame()

    log.info("Control conditions: %s", ctrl)
    log.info("Disease conditions: %s", dis)

    # Top 5 cell types by cell count
    ct_counts = adata.obs[ct_col].value_counts()
    top_cts = ct_counts.head(5).index.tolist()
    log.info("Running DE for top %d cell types: %s", len(top_cts), top_cts)

    all_de = []
    for ct in top_cts:
        ct_mask = (adata.obs[ct_col] == ct).values
        obs_ct = adata.obs.loc[ct_mask, [sample_col, cond_col]].copy()
        X_ct = X_lncrna[ct_mask, :]

        # Aggregate to pseudobulk (sum per sample)
        obs_ct = obs_ct.reset_index(drop=True)
        samples = obs_ct[sample_col].unique()

        pb_data = []
        pb_cond = []
        for samp in samples:
            samp_mask = (obs_ct[sample_col] == samp).values
            n_samp_cells = samp_mask.sum()
            if n_samp_cells < 3:
                continue
            samp_sum = X_ct[samp_mask, :].sum(axis=0)
            pb_data.append(np.asarray(samp_sum).ravel())
            # Condition for this sample (take most common)
            conds = obs_ct.loc[samp_mask, cond_col]
            pb_cond.append(conds.mode().iloc[0] if len(conds) > 0 else "Unknown")

        if len(pb_data) < 4:
            log.info("  %s: insufficient pseudobulk samples (%d) — skipping", ct, len(pb_data))
            continue

        pb_mat = np.vstack(pb_data)  # (n_samples, n_genes)
        pb_cond = np.array(pb_cond)

        is_ctrl = np.isin(pb_cond, ctrl)
        is_dis = np.isin(pb_cond, dis)

        n_ctrl = is_ctrl.sum()
        n_dis = is_dis.sum()
        log.info("  %s: %d ctrl, %d disease pseudobulk samples", ct, n_ctrl, n_dis)

        if n_ctrl < 2 or n_dis < 2:
            log.info("  %s: insufficient samples per group — skipping", ct)
            continue

        # Filter genes: detected in >= 5% of pseudobulk samples
        det = (pb_mat > 0).mean(axis=0)
        gene_mask = det >= 0.05
        n_tested = gene_mask.sum()
        log.info("  %s: testing %d / %d lncRNA genes (>= 5%% detection)", ct, n_tested, len(lncrna_genes))

        pvals = []
        logfcs = []
        mean_dis_list = []
        mean_ctrl_list = []
        tested_genes = []

        for g_idx in range(len(lncrna_genes)):
            if not gene_mask[g_idx]:
                continue
            x_ctrl = pb_mat[is_ctrl, g_idx]
            x_dis = pb_mat[is_dis, g_idx]

            # Log2 fold change (pseudocount 1)
            mean_c = x_ctrl.mean()
            mean_d = x_dis.mean()
            lfc = np.log2((mean_d + 1) / (mean_c + 1))

            try:
                stat, pval = mannwhitneyu(x_dis, x_ctrl, alternative="two-sided")
            except ValueError:
                pval = 1.0

            pvals.append(pval)
            logfcs.append(lfc)
            mean_dis_list.append(mean_d)
            mean_ctrl_list.append(mean_c)
            tested_genes.append(lncrna_genes[g_idx])

        if len(pvals) == 0:
            continue

        # BH correction
        _, padj, _, _ = multipletests(pvals, method="fdr_bh")

        de_df = pd.DataFrame({
            "gene": tested_genes,
            "cell_type": ct,
            "mean_disease": np.round(mean_dis_list, 4),
            "mean_control": np.round(mean_ctrl_list, 4),
            "logFC": np.round(logfcs, 4),
            "pval": pvals,
            "padj": padj,
        })
        all_de.append(de_df)

        n_sig = (de_df["padj"] < 0.1).sum()
        log.info("  %s: %d DE lncRNAs at padj < 0.1", ct, n_sig)

    if all_de:
        result = pd.concat(all_de, ignore_index=True)
    else:
        result = pd.DataFrame(columns=[
            "gene", "cell_type", "mean_disease", "mean_control",
            "logFC", "pval", "padj",
        ])

    elapsed = time.time() - t0
    log.info("Pseudobulk DE complete in %.1f s — %d total rows", elapsed, len(result))
    return result


# ============================================================================
# Section 7 — Cross-reference with SCENIC+ regulons
# ============================================================================
def scenic_cross_reference(tau_df, scenic_path):
    """Check overlap of cell-type-specific lncRNAs with SCENIC+ regulon targets."""
    log.info("=" * 70)
    log.info("SECTION 7: Cross-reference with SCENIC+ regulons")
    log.info("=" * 70)

    specific = set(tau_df.loc[tau_df["tau"] > 0.8, "gene"]) if len(tau_df) > 0 else set()
    log.info("Cell-type-specific lncRNAs (tau > 0.8): %d", len(specific))

    if not os.path.exists(scenic_path):
        log.warning("SCENIC+ regulon file not found: %s", scenic_path)
        return 0, set()

    try:
        reg = pd.read_csv(scenic_path)
        regulon_targets = set(reg["target_gene"].dropna().unique())
        log.info("SCENIC+ regulon targets: %d unique genes", len(regulon_targets))

        overlap = specific & regulon_targets
        log.info("Overlap (specific lncRNAs in SCENIC+ regulon targets): %d", len(overlap))
        if overlap:
            log.info("Overlapping genes: %s", sorted(overlap)[:20])
        return len(overlap), overlap
    except Exception as e:
        log.warning("Failed to load SCENIC+ regulons: %s", e)
        return 0, set()


# ============================================================================
# Section 8 — Summary and output
# ============================================================================
def write_outputs(out_dir, expr_df, tau_df, de_df, qc_df):
    """Write all output CSVs."""
    log.info("=" * 70)
    log.info("SECTION 8: Writing outputs to %s", out_dir)
    log.info("=" * 70)
    os.makedirs(out_dir, exist_ok=True)

    expr_path = os.path.join(out_dir, "lncrna_celltype_expression.csv")
    expr_df.to_csv(expr_path, index=False)
    log.info("  %s — %d rows", expr_path, len(expr_df))

    tau_path = os.path.join(out_dir, "lncrna_celltype_specific.csv")
    tau_df.to_csv(tau_path, index=False)
    log.info("  %s — %d rows", tau_path, len(tau_df))

    de_path = os.path.join(out_dir, "lncrna_celltype_de.csv")
    de_df.to_csv(de_path, index=False)
    log.info("  %s — %d rows", de_path, len(de_df))

    qc_path = os.path.join(out_dir, "elatus_qc_metrics.csv")
    qc_df.to_csv(qc_path, index=False)
    log.info("  %s — %d rows", qc_path, len(qc_df))


def print_summary(expr_df, tau_df, de_df, qc_df, n_scenic_overlap):
    """Print final summary statistics."""
    log.info("=" * 70)
    log.info("SUMMARY")
    log.info("=" * 70)

    n_genes = tau_df["gene"].nunique() if len(tau_df) > 0 else 0
    n_specific = (tau_df["tau"] > 0.8).sum() if len(tau_df) > 0 else 0
    n_ubiq = (tau_df["tau"] < 0.3).sum() if len(tau_df) > 0 else 0
    n_de_sig = (de_df["padj"] < 0.1).sum() if len(de_df) > 0 else 0
    n_de_strict = (de_df["padj"] < 0.05).sum() if len(de_df) > 0 else 0
    n_celltypes_de = de_df["cell_type"].nunique() if len(de_df) > 0 else 0

    mean_cf = qc_df["correction_factor"].mean() if len(qc_df) > 0 else 0
    median_frac = qc_df["lncrna_frac_observed"].median() if len(qc_df) > 0 and "lncrna_frac_observed" in qc_df.columns else 0

    log.info("lncRNA genes profiled:            %d", n_genes)
    log.info("Cell-type-specific (tau > 0.8):   %d", n_specific)
    log.info("Ubiquitous (tau < 0.3):           %d", n_ubiq)
    log.info("DE lncRNAs (padj < 0.1):          %d across %d cell types", n_de_sig, n_celltypes_de)
    log.info("DE lncRNAs (padj < 0.05):         %d", n_de_strict)
    log.info("SCENIC+ regulon overlap:          %d", n_scenic_overlap)
    log.info("Mean correction factor:           %.2f (1.0 = no correction)", mean_cf)
    log.info("Median lncRNA fraction:           %.4f", median_frac)
    log.info("=" * 70)


# ============================================================================
# Main
# ============================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Cell-type-specific ncRNA profiling (ELATUS-inspired correction)",
    )
    parser.add_argument(
        "--atlas", default=DEFAULT_ATLAS,
        help="Path to integrated scRNA-seq atlas h5ad",
    )
    parser.add_argument(
        "--output-dir", default=DEFAULT_OUT,
        help="Output directory for results",
    )
    parser.add_argument(
        "--gtf", default=GENCODE_GTF,
        help="Path to GENCODE GTF (gzipped)",
    )
    parser.add_argument(
        "--dream", default=DREAM_RESULTS,
        help="Path to dream bulk DE results CSV",
    )
    parser.add_argument(
        "--scenic", default=SCENIC_REGULONS,
        help="Path to SCENIC+ hepatocyte regulons CSV",
    )
    args = parser.parse_args()

    start_time = datetime.now()
    log.info("=" * 70)
    log.info("ELATUS-inspired lncRNA profiling — %s", start_time.strftime("%Y-%m-%d %H:%M:%S"))
    log.info("=" * 70)
    log.info("Atlas:      %s", args.atlas)
    log.info("GTF:        %s", args.gtf)
    log.info("Dream:      %s", args.dream)
    log.info("Output:     %s", args.output_dir)

    # Section 1: Parse GTF
    gene_biotype, biotype_sets = parse_gencode_biotypes(args.gtf)
    lncrna_set = biotype_sets["lncRNA"]

    # Section 2: Load atlas
    adata, ct_col, cond_col, sample_col = load_atlas(args.atlas)

    if not ct_col:
        log.error("Could not identify cell type column — aborting")
        sys.exit(1)

    # Identify lncRNA genes in atlas
    lncrna_mask, lncrna_matched, symbol_col = identify_lncrna_genes(adata, lncrna_set)
    n_lncrna = lncrna_mask.sum()
    log.info("lncRNA genes in atlas: %d", n_lncrna)

    if n_lncrna == 0:
        log.error("No lncRNA genes found in atlas — check gene naming. Aborting.")
        sys.exit(1)

    # Get lncRNA gene names in atlas order
    if symbol_col:
        lncrna_genes = adata.var.loc[lncrna_mask, symbol_col].tolist()
    else:
        lncrna_genes = adata.var_names[lncrna_mask].tolist()

    # Section 3: ELATUS correction
    X_lncrna, qc_df, correction_factors = elatus_correction(
        adata, lncrna_mask, ct_col, args.dream
    )

    # Section 4: Expression profiling
    expr_df = celltype_expression(adata, X_lncrna, lncrna_genes, ct_col, cond_col)

    # Section 5: Tau specificity
    tau_df = compute_tau(adata, X_lncrna, lncrna_genes, ct_col)

    # Section 6: Pseudobulk DE
    de_df = pseudobulk_de(adata, X_lncrna, lncrna_genes, ct_col, cond_col, sample_col)

    # Section 7: SCENIC+ cross-reference
    n_scenic_overlap, scenic_genes = scenic_cross_reference(tau_df, args.scenic)

    # Section 8: Write outputs and summary
    write_outputs(args.output_dir, expr_df, tau_df, de_df, qc_df)
    print_summary(expr_df, tau_df, de_df, qc_df, n_scenic_overlap)

    elapsed = datetime.now() - start_time
    log.info("Total runtime: %s", str(elapsed).split(".")[0])
    log.info("Done.")


if __name__ == "__main__":
    main()
