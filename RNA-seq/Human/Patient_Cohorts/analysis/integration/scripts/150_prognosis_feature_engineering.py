#!/usr/bin/env python3
"""
150_prognosis_feature_engineering.py — Prognosis Feature Engineering
--------------------------------------------------------------------
Assembles all features and pseudo-labels for the Progression Risk Score
pipeline. Produces a ~255-column feature matrix (1,444 samples) with 7
feature tiers and multiple pseudo-label variants.

Feature tiers:
  T1: Divergence gene expression (top 100 by |cohens_d|) — 100 features
  T2: Cell-type fractions (13 BayesPrism + 4 engineered) — 17 features
  T3: Latent embeddings (NAS-VAE 64-d + BulkFormer 64-d) — 128 features
  T4: GAS6::MERTK interaction features — 2 features
  T5: Transition program scores (4 fibrosis transitions) — 4 features
  T6: HCC genetic risk score — 1 feature
  T7: Clinical covariates (sex, fibrosis_stage, age) — 3 features

Pseudo-labels:
  - CPS (Composite Progression Score): weighted combination
  - Binary S1/S2 (NMF subtypes)
  - Continuous P(F4) (fate probability)
  - Residualized P(F4) (progression beyond stage)

Input files:
  - results/progression/divergence_genes.csv
  - results/staging_classifier/prepared_data.h5  (rank_expression 3000x1444)
  - results/progression/cibersortx_celltype_expression/bayesprism_proportions.csv
  - results/staging_classifier/nas_embeddings_all_samples.csv
  - results/staging_classifier/bulkformer_embeddings.csv
  - results/progression/consensus_pseudotime.csv
  - results/progression/fate_probabilities.csv
  - results/progression/differentiation_potential.csv
  - results/staging_classifier/modeling_metadata.csv
  - results/progression/driver_scores.csv
  - RNA-seq/results/subtypes/nmf_assignments.csv
  - RNA-seq/results/causal_inference/finngen_hcc/coloc_results.csv
  - RNA-seq/results/causal_inference/ghouse_hcc/coloc_results.csv
  - results/gene_annotation/human_ensg_to_symbol.tsv

Output (to results/prognosis/):
  - prognosis_feature_matrix.csv      (1,444 x ~255 cols)
  - prognosis_pseudo_labels.csv       (1,444 x label cols)
  - prognosis_feature_description.csv (feature_name, tier, description)
  - hcc_coloc_genes.csv               (genes used for HCC score)

SLURM: io, 4 CPUs, 16G RAM, 48h
Env:   micromamba activate spatial
"""

import logging
import os
import sys
import time
import warnings

import h5py
import numpy as np
import pandas as pd
from scipy.spatial import KDTree

warnings.filterwarnings("ignore", category=RuntimeWarning)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INTEG = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RNASEQ = os.path.join(BASE, "RNA-seq")

# Input files
DIVERGENCE_GENES = os.path.join(INTEG, "results/progression/divergence_genes.csv")
PREPARED_H5 = os.path.join(INTEG, "results/staging_classifier/prepared_data.h5")
BAYESPRISM_PROP = os.path.join(
    INTEG, "results/progression/cibersortx_celltype_expression/bayesprism_proportions.csv"
)
NAS_EMBEDDINGS = os.path.join(INTEG, "results/staging_classifier/nas_embeddings_all_samples.csv")
BULKFORMER_EMBEDDINGS = os.path.join(INTEG, "results/staging_classifier/bulkformer_embeddings.csv")
CONSENSUS_PT = os.path.join(INTEG, "results/progression/consensus_pseudotime.csv")
FATE_PROBS = os.path.join(INTEG, "results/progression/fate_probabilities.csv")
DIFF_POTENTIAL = os.path.join(INTEG, "results/progression/differentiation_potential.csv")
MODELING_META = os.path.join(INTEG, "results/staging_classifier/modeling_metadata.csv")
DRIVER_SCORES = os.path.join(INTEG, "results/progression/driver_scores.csv")
NMF_ASSIGNMENTS = os.path.join(RNASEQ, "results/subtypes/nmf_assignments.csv")
GENE_ANNOT = os.path.join(INTEG, "results/gene_annotation/human_ensg_to_symbol.tsv")
FINNGEN_HCC_COLOC = os.path.join(RNASEQ, "results/causal_inference/finngen_hcc/coloc_results.csv")
GHOUSE_HCC_COLOC = os.path.join(RNASEQ, "results/causal_inference/ghouse_hcc/coloc_results.csv")

# Output directory
OUT_DIR = os.path.join(INTEG, "results/prognosis")

# ── Parameters ───────────────────────────────────────────────────────────────
N_DIVERGENCE_GENES = 100  # top genes by |cohens_d|
N_DRIVER_GENES_PER_TRANSITION = 50  # top driver genes per transition
HCC_COLOC_THRESHOLD = 0.3  # PP.H4 threshold for HCC genetic score
KNN_K = 30  # k for pseudotime velocity (nearest neighbors)
FIBROSIS_TRANSITIONS = ["F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4"]

# CPS weights
CPS_W_S2 = 0.25
CPS_W_PF4 = 0.30
CPS_W_DIFFPOT = 0.20
CPS_W_STELLATE = 0.15
CPS_W_VELOCITY = 0.10


# ── Utility ──────────────────────────────────────────────────────────────────
def safe_load(path, reader="csv", **kwargs):
    """Load a file with existence check. Returns None if missing."""
    if not os.path.exists(path):
        log.warning("File not found: %s", path)
        return None
    log.info("Loading %s", os.path.basename(path))
    if reader == "csv":
        return pd.read_csv(path, **kwargs)
    elif reader == "tsv":
        return pd.read_csv(path, sep="\t", **kwargs)
    return None


def minmax_scale(x):
    """Min-max scale to [0, 1], handling NaN and all-NaN arrays."""
    x = np.asarray(x, dtype=float)
    if np.all(np.isnan(x)):
        return x.copy()  # all NaN -> return NaN
    xmin = np.nanmin(x)
    xmax = np.nanmax(x)
    denom = xmax - xmin
    if denom == 0:
        return np.where(np.isnan(x), np.nan, 0.5)
    result = (x - xmin) / denom
    return result


# ── Main ─────────────────────────────────────────────────────────────────────
def main():
    t0 = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    # ------------------------------------------------------------------
    # 0. Load master sample list and gene annotation
    # ------------------------------------------------------------------
    meta = safe_load(MODELING_META)
    if meta is None:
        log.error("modeling_metadata.csv is required. Exiting.")
        sys.exit(1)
    sample_ids = meta["sample_id"].values
    n_samples = len(sample_ids)
    log.info("Master sample list: %d samples", n_samples)

    gene_annot = safe_load(GENE_ANNOT, reader="tsv")
    if gene_annot is None:
        log.error("Gene annotation file is required. Exiting.")
        sys.exit(1)
    # Build mapping: gene_base (no version) -> symbol, and versioned -> symbol
    ensg_to_symbol = dict(zip(gene_annot["gene_base"], gene_annot["symbol"]))
    ensg_versioned_to_symbol = dict(zip(gene_annot["gene_id"], gene_annot["symbol"]))
    symbol_to_ensg_versioned = {}
    for gid, sym in ensg_versioned_to_symbol.items():
        if sym not in symbol_to_ensg_versioned:
            symbol_to_ensg_versioned[sym] = gid
    log.info("Gene annotation: %d versioned, %d base mappings", len(ensg_versioned_to_symbol), len(ensg_to_symbol))

    # ------------------------------------------------------------------
    # 1. Load h5 expression data
    # ------------------------------------------------------------------
    log.info("Loading prepared_data.h5 ...")
    h5 = h5py.File(PREPARED_H5, "r")
    h5_gene_names = np.array([g.decode() for g in h5["gene_names"][:]])  # versioned Ensembl
    h5_sample_ids = np.array([s.decode() for s in h5["sample_ids"][:]])
    rank_expr = h5["rank_expression"][:].T  # (1444, 3000) after transpose
    h5.close()
    log.info("h5 rank expression: %s (samples x genes)", rank_expr.shape)

    # Build gene lookup: versioned Ensembl -> column index in h5
    h5_gene_idx = {g: i for i, g in enumerate(h5_gene_names)}
    # Also build symbol -> column index
    h5_symbol_idx = {}
    for versioned, idx in h5_gene_idx.items():
        sym = ensg_versioned_to_symbol.get(versioned)
        if sym:
            h5_symbol_idx[sym] = idx

    # Verify sample order matches metadata
    assert np.array_equal(h5_sample_ids, sample_ids), (
        "Sample order mismatch between h5 and modeling_metadata.csv"
    )

    # Feature matrix and descriptions accumulator
    feature_dict = {}  # feature_name -> np.array(n_samples)
    desc_rows = []  # (feature_name, tier, description)

    # ==================================================================
    # TIER 1: Divergence gene expression (100 features)
    # ==================================================================
    log.info("=== TIER 1: Divergence gene expression ===")
    div_genes = safe_load(DIVERGENCE_GENES)
    t1_count = 0
    if div_genes is not None:
        div_genes["abs_cohens_d"] = div_genes["cohens_d"].abs()
        top_div = div_genes.nlargest(N_DIVERGENCE_GENES, "abs_cohens_d")
        for _, row in top_div.iterrows():
            sym = row["gene_symbol"]
            ensg = row["gene"]  # versioned Ensembl ID
            # Try versioned ID first, then symbol lookup
            col_idx = h5_gene_idx.get(ensg)
            if col_idx is None:
                col_idx = h5_symbol_idx.get(sym)
            if col_idx is not None:
                fname = f"div_{sym}"
                feature_dict[fname] = rank_expr[:, col_idx].copy()
                desc_rows.append((fname, "T1_divergence", f"Rank expression of {sym} (|d|={row['abs_cohens_d']:.3f})"))
                t1_count += 1
            else:
                log.warning("T1: gene %s (%s) not in h5 expression matrix", sym, ensg)
        log.info("T1: %d / %d divergence genes mapped", t1_count, N_DIVERGENCE_GENES)
    else:
        log.warning("T1: divergence_genes.csv not found, skipping")

    # ==================================================================
    # TIER 2: Cell-type fractions (17 features)
    # ==================================================================
    log.info("=== TIER 2: Cell-type fractions ===")
    bp_prop = safe_load(BAYESPRISM_PROP)
    t2_count = 0
    if bp_prop is not None:
        bp_prop = bp_prop.set_index("sample_id").reindex(sample_ids)
        celltype_cols = [c for c in bp_prop.columns if c != "sample_id"]
        log.info("T2: %d cell types from BayesPrism", len(celltype_cols))

        # 13 raw proportions
        for ct in celltype_cols:
            fname = f"ct_{ct}"
            feature_dict[fname] = bp_prop[ct].values.astype(float)
            desc_rows.append((fname, "T2_celltype", f"BayesPrism fraction: {ct}"))
            t2_count += 1

        # Engineered features
        stellate = bp_prop.get("Stellate", pd.Series(np.nan, index=bp_prop.index))
        hepatocyte = bp_prop.get("Hepatocyte", pd.Series(np.nan, index=bp_prop.index))
        macrophage = bp_prop.get("Macrophage", pd.Series(np.nan, index=bp_prop.index))
        monocyte = bp_prop.get("Monocyte", pd.Series(np.nan, index=bp_prop.index))
        nk = bp_prop.get("NK_cell", pd.Series(np.nan, index=bp_prop.index))
        tcell = bp_prop.get("T_cell", pd.Series(np.nan, index=bp_prop.index))
        cholangio = bp_prop.get("Cholangiocyte", pd.Series(np.nan, index=bp_prop.index))

        # stellate_gt4pct (binary)
        stellate_binary = (stellate > 0.04).astype(float)
        stellate_binary[stellate.isna()] = np.nan
        feature_dict["ct_stellate_gt4pct"] = stellate_binary.values
        desc_rows.append(("ct_stellate_gt4pct", "T2_celltype", "Binary: stellate fraction > 4%"))
        t2_count += 1

        # stellate_hepatocyte_ratio
        ratio = stellate / (hepatocyte + 1e-6)
        feature_dict["ct_stellate_hepatocyte_ratio"] = ratio.values.astype(float)
        desc_rows.append(("ct_stellate_hepatocyte_ratio", "T2_celltype", "Stellate / (Hepatocyte + 1e-6) ratio"))
        t2_count += 1

        # immune_infiltration
        immune = macrophage + monocyte + nk + tcell
        feature_dict["ct_immune_infiltration"] = immune.values.astype(float)
        desc_rows.append(("ct_immune_infiltration", "T2_celltype", "Sum of Macrophage + Monocyte + NK_cell + T_cell"))
        t2_count += 1

        # fibrogenic_axis
        fibro = stellate + cholangio
        feature_dict["ct_fibrogenic_axis"] = fibro.values.astype(float)
        desc_rows.append(("ct_fibrogenic_axis", "T2_celltype", "Sum of Stellate + Cholangiocyte"))
        t2_count += 1

        log.info("T2: %d features total", t2_count)
    else:
        log.warning("T2: bayesprism_proportions.csv not found, skipping")

    # ==================================================================
    # TIER 3: Embeddings (128 features)
    # ==================================================================
    log.info("=== TIER 3: Latent embeddings ===")
    t3_count = 0

    # NAS-VAE (64-dim)
    nas_emb = safe_load(NAS_EMBEDDINGS)
    if nas_emb is not None:
        nas_emb = nas_emb.set_index("sample_id").reindex(sample_ids)
        z_cols = [c for c in nas_emb.columns if c.startswith("z")]
        for c in z_cols:
            fname = f"nasvae_{c}"
            feature_dict[fname] = nas_emb[c].values.astype(float)
            desc_rows.append((fname, "T3_embedding", f"NAS-VAE latent dim {c}"))
            t3_count += 1
        log.info("T3: %d NAS-VAE dimensions", len(z_cols))
    else:
        log.warning("T3: NAS-VAE embeddings not found, skipping")

    # BulkFormer (64-dim)
    bf_emb = safe_load(BULKFORMER_EMBEDDINGS)
    if bf_emb is not None:
        bf_emb = bf_emb.set_index("sample_id").reindex(sample_ids)
        bf_cols = [c for c in bf_emb.columns if c.startswith("bf")]
        for c in bf_cols:
            fname = f"bulkformer_{c}"
            feature_dict[fname] = bf_emb[c].values.astype(float)
            desc_rows.append((fname, "T3_embedding", f"BulkFormer latent dim {c}"))
            t3_count += 1
        log.info("T3: %d BulkFormer dimensions", len(bf_cols))
    else:
        log.warning("T3: BulkFormer embeddings not found, skipping")

    log.info("T3: %d embedding features total", t3_count)

    # ==================================================================
    # TIER 4: GAS6::MERTK interaction features (2 features)
    # ==================================================================
    log.info("=== TIER 4: GAS6::MERTK features ===")
    t4_count = 0
    gas6_idx = h5_symbol_idx.get("GAS6")
    mertk_idx = h5_symbol_idx.get("MERTK")

    if gas6_idx is not None and mertk_idx is not None:
        gas6_expr = rank_expr[:, gas6_idx]
        mertk_expr = rank_expr[:, mertk_idx]

        # Product
        feature_dict["gas6_mertk_product"] = gas6_expr * mertk_expr
        desc_rows.append(("gas6_mertk_product", "T4_gas6_mertk", "GAS6 * MERTK rank expression product"))
        t4_count += 1

        # Ratio
        feature_dict["gas6_mertk_ratio"] = gas6_expr / (mertk_expr + 1e-6)
        desc_rows.append(("gas6_mertk_ratio", "T4_gas6_mertk", "GAS6 / (MERTK + 1e-6) rank expression ratio"))
        t4_count += 1

        log.info("T4: %d GAS6::MERTK features", t4_count)
    else:
        missing = []
        if gas6_idx is None:
            missing.append("GAS6")
        if mertk_idx is None:
            missing.append("MERTK")
        log.warning("T4: %s not found in h5, skipping", ", ".join(missing))

    # ==================================================================
    # TIER 5: Transition program scores (4 features)
    # ==================================================================
    log.info("=== TIER 5: Transition program scores ===")
    driver_df = safe_load(DRIVER_SCORES, low_memory=False)
    t5_count = 0
    if driver_df is not None:
        for trans in FIBROSIS_TRANSITIONS:
            sub = driver_df[driver_df["transition"] == trans].copy()
            if len(sub) == 0:
                log.warning("T5: no genes for transition %s", trans)
                continue
            # Top N driver genes by driver_score
            top_drivers = sub.nlargest(N_DRIVER_GENES_PER_TRANSITION, "driver_score")
            # Map gene symbols to h5 indices
            gene_indices = []
            for sym in top_drivers["gene_symbol"]:
                idx = h5_symbol_idx.get(sym)
                if idx is not None:
                    gene_indices.append(idx)
            if len(gene_indices) == 0:
                log.warning("T5: no driver genes mapped for %s", trans)
                continue
            # Mean rank expression across top driver genes
            scores = rank_expr[:, gene_indices].mean(axis=1)
            fname = f"trans_{trans}"
            feature_dict[fname] = scores
            desc_rows.append((
                fname, "T5_transition",
                f"Mean rank expr of top {len(gene_indices)} drivers for {trans}"
            ))
            t5_count += 1
            log.info("T5: %s — %d / %d driver genes mapped", trans, len(gene_indices), N_DRIVER_GENES_PER_TRANSITION)
        log.info("T5: %d transition features", t5_count)
    else:
        log.warning("T5: driver_scores.csv not found, skipping")

    # ==================================================================
    # TIER 6: HCC genetic risk score (1 feature)
    # ==================================================================
    log.info("=== TIER 6: HCC genetic risk score ===")
    hcc_coloc_all = []
    for path, source in [(FINNGEN_HCC_COLOC, "finngen_hcc"), (GHOUSE_HCC_COLOC, "ghouse_hcc")]:
        df = safe_load(path)
        if df is not None:
            df = df[["gene", "PP.H4"]].copy()
            df["source"] = source
            hcc_coloc_all.append(df)
            log.info("T6: %d genes from %s", len(df), source)

    t6_count = 0
    hcc_genes_out = pd.DataFrame()
    if len(hcc_coloc_all) > 0:
        hcc_combined = pd.concat(hcc_coloc_all, ignore_index=True)
        # Take max PP.H4 per gene across sources
        hcc_max = hcc_combined.groupby("gene")["PP.H4"].max().reset_index()
        hcc_sig = hcc_max[hcc_max["PP.H4"] > HCC_COLOC_THRESHOLD].copy()
        log.info("T6: %d genes with PP.H4 > %.2f", len(hcc_sig), HCC_COLOC_THRESHOLD)

        if len(hcc_sig) > 0:
            # Map genes to h5 indices and compute weighted z-scored expression
            hcc_gene_weights = []
            hcc_gene_indices = []
            for _, row in hcc_sig.iterrows():
                sym = row["gene"]
                idx = h5_symbol_idx.get(sym)
                if idx is not None:
                    hcc_gene_indices.append(idx)
                    hcc_gene_weights.append(row["PP.H4"])
            if len(hcc_gene_indices) > 0:
                weights = np.array(hcc_gene_weights)
                # Z-score each gene across samples, then weighted sum
                gene_exprs = rank_expr[:, hcc_gene_indices]  # (n_samples, n_genes)
                gene_z = np.zeros_like(gene_exprs)
                for j in range(gene_exprs.shape[1]):
                    col = gene_exprs[:, j]
                    mu = np.nanmean(col)
                    sd = np.nanstd(col)
                    if sd > 0:
                        gene_z[:, j] = (col - mu) / sd
                    else:
                        gene_z[:, j] = 0.0
                hcc_score = gene_z @ weights
                feature_dict["hcc_genetic_score"] = hcc_score
                desc_rows.append((
                    "hcc_genetic_score", "T6_hcc",
                    f"Weighted sum of z-scored expr for {len(hcc_gene_indices)} HCC COLOC genes (PP.H4 weights)"
                ))
                t6_count = 1
                log.info("T6: HCC score from %d / %d genes", len(hcc_gene_indices), len(hcc_sig))

                # Save HCC gene details
                hcc_genes_out = hcc_sig.copy()
                hcc_genes_out["in_expression_matrix"] = hcc_genes_out["gene"].isin(
                    [ensg_versioned_to_symbol.get(h5_gene_names[i], h5_gene_names[i])
                     for i in hcc_gene_indices]
                )
            else:
                log.warning("T6: no HCC COLOC genes mapped to expression matrix")
        else:
            log.warning("T6: no genes pass PP.H4 threshold")
    else:
        log.warning("T6: no HCC COLOC files found, skipping")

    # ==================================================================
    # TIER 7: Clinical covariates (3 features)
    # ==================================================================
    log.info("=== TIER 7: Clinical covariates ===")
    # sex: M=0, F=1, missing=0.5
    sex_map = {"M": 0.0, "F": 1.0}
    sex_vals = meta["sex"].map(sex_map).fillna(0.5).values.astype(float)
    feature_dict["clin_sex"] = sex_vals
    desc_rows.append(("clin_sex", "T7_clinical", "Sex (M=0, F=1, missing=0.5)"))

    # fibrosis_stage: ordinal 0-4, -1 for missing
    fib_vals = meta["fibrosis_stage"].copy()
    fib_vals = pd.to_numeric(fib_vals, errors="coerce")
    fib_vals = fib_vals.fillna(-1).values.astype(float)
    feature_dict["clin_fibrosis_stage"] = fib_vals
    desc_rows.append(("clin_fibrosis_stage", "T7_clinical", "Fibrosis stage (ordinal 0-4, -1=missing)"))

    # age: continuous, NaN for missing
    age_vals = meta["age"].copy()
    age_vals = pd.to_numeric(age_vals, errors="coerce").values.astype(float)
    feature_dict["clin_age"] = age_vals
    desc_rows.append(("clin_age", "T7_clinical", "Age (continuous, NaN=missing)"))

    log.info("T7: 3 clinical covariates")

    # ------------------------------------------------------------------
    # Assemble feature matrix
    # ------------------------------------------------------------------
    log.info("Assembling feature matrix ...")
    feat_df = pd.DataFrame(feature_dict, index=sample_ids)
    feat_df.index.name = "sample_id"
    n_features = feat_df.shape[1]
    log.info("Feature matrix: %d samples x %d features", n_samples, n_features)

    # Tier summary
    tier_counts = {}
    for _, tier, _ in desc_rows:
        tier_counts[tier] = tier_counts.get(tier, 0) + 1
    for tier, cnt in sorted(tier_counts.items()):
        log.info("  %s: %d features", tier, cnt)

    # ==================================================================
    # Build pseudo-labels
    # ==================================================================
    log.info("=== Building pseudo-labels ===")
    labels = pd.DataFrame(index=sample_ids)
    labels.index.name = "sample_id"

    # Load auxiliary data (index on sample_id, reindex to master list)
    nmf = safe_load(NMF_ASSIGNMENTS)
    fate = safe_load(FATE_PROBS)
    cpt = safe_load(CONSENSUS_PT)
    dpot = safe_load(DIFF_POTENTIAL)

    # S2 binary (1,254 samples; rest NaN)
    if nmf is not None:
        nmf_map = nmf.set_index("sample_id")["nmf_subtype"].reindex(sample_ids)
        labels["s2_binary"] = (nmf_map == "S2").astype(float)
        labels.loc[nmf_map.isna(), "s2_binary"] = np.nan
        log.info("Pseudo-label: s2_binary — %d non-NaN", labels["s2_binary"].notna().sum())
    else:
        labels["s2_binary"] = np.nan
        log.warning("NMF assignments not found, s2_binary all NaN")

    # P(F4) (1,254 samples; rest NaN)
    if fate is not None:
        fate_map = fate.set_index("sample_id")["fate_prob_F4"].reindex(sample_ids)
        labels["p_f4"] = fate_map.values.astype(float)
        log.info("Pseudo-label: p_f4 — %d non-NaN", labels["p_f4"].notna().sum())
    else:
        labels["p_f4"] = np.nan
        log.warning("Fate probabilities not found, p_f4 all NaN")

    # Residualized P(F4): regress out fibrosis_stage
    if fate is not None:
        pf4 = labels["p_f4"].values.copy()
        fib_for_resid = fib_vals.copy()
        # Only use samples with valid P(F4) and valid fibrosis stage
        valid = (~np.isnan(pf4)) & (fib_for_resid >= 0)
        if valid.sum() > 10:
            x_valid = fib_for_resid[valid]
            y_valid = pf4[valid]
            # Linear regression: P(F4) ~ fibrosis_stage
            coeffs = np.polyfit(x_valid, y_valid, deg=1)
            predicted = np.polyval(coeffs, fib_for_resid)
            residual = pf4 - predicted
            residual[~valid] = np.nan
            labels["p_f4_residualized"] = residual
            log.info("Pseudo-label: p_f4_residualized — %d non-NaN (slope=%.4f)",
                      (~np.isnan(residual)).sum(), coeffs[0])
        else:
            labels["p_f4_residualized"] = np.nan
            log.warning("Too few valid samples for P(F4) residualization")
    else:
        labels["p_f4_residualized"] = np.nan

    # Pseudotime and differentiation potential (1,444 samples)
    pseudotime_vals = np.full(n_samples, np.nan)
    if cpt is not None:
        pt_map = cpt.set_index("sample_id")["pseudotime_consensus"].reindex(sample_ids)
        pseudotime_vals = pt_map.values.astype(float)
        log.info("Pseudotime: %d non-NaN", (~np.isnan(pseudotime_vals)).sum())

    diffpot_vals = np.full(n_samples, np.nan)
    if dpot is not None:
        dp_map = dpot.set_index("sample_id")["differentiation_potential"].reindex(sample_ids)
        diffpot_vals = dp_map.values.astype(float)
        log.info("Diff potential: %d non-NaN", (~np.isnan(diffpot_vals)).sum())

    # Pseudotime velocity: |PT(sample) - mean(PT(k nearest neighbors in NAS-VAE space))|
    velocity_vals = np.full(n_samples, np.nan)
    if nas_emb is not None and cpt is not None:
        log.info("Computing pseudotime velocity (k=%d NAS-VAE neighbors) ...", KNN_K)
        # nas_emb was already loaded and reindexed to sample_ids in Tier 3
        z_cols_vel = [c for c in nas_emb.columns if c.startswith("z")]
        emb_matrix = nas_emb[z_cols_vel].values.astype(float)

        # Only compute for samples with valid embeddings and pseudotime
        valid_emb = ~np.isnan(emb_matrix).any(axis=1)
        valid_pt = ~np.isnan(pseudotime_vals)
        valid_both = valid_emb & valid_pt
        n_valid = valid_both.sum()
        log.info("Velocity: %d samples with valid embeddings + pseudotime", n_valid)

        if n_valid > KNN_K + 1:
            # Build KDTree on valid samples
            valid_indices = np.where(valid_both)[0]
            tree = KDTree(emb_matrix[valid_both])
            pt_valid = pseudotime_vals[valid_both]

            # Query k+1 neighbors (first is self)
            dists, nn_idx = tree.query(emb_matrix[valid_both], k=KNN_K + 1)
            # Mean pseudotime of k nearest neighbors (exclude self = index 0)
            neighbor_pt = pt_valid[nn_idx[:, 1:]]  # (n_valid, k)
            mean_neighbor_pt = np.nanmean(neighbor_pt, axis=1)
            local_velocity = np.abs(pt_valid - mean_neighbor_pt)

            # Map back to full array
            for i, full_idx in enumerate(valid_indices):
                velocity_vals[full_idx] = local_velocity[i]
            log.info("Velocity: computed for %d samples (mean=%.4f, std=%.4f)",
                      n_valid, np.nanmean(velocity_vals), np.nanstd(velocity_vals))
        else:
            log.warning("Velocity: too few valid samples (%d <= k=%d)", n_valid, KNN_K)

    # Composite Progression Score (CPS)
    log.info("Computing CPS ...")
    stellate_vals = feature_dict.get("ct_Stellate")
    if stellate_vals is None:
        stellate_vals = np.full(n_samples, np.nan)

    # Scale each component to [0,1]
    s2_scaled = labels["s2_binary"].values.astype(float)  # already 0/1
    pf4_scaled = minmax_scale(labels["p_f4"].values.astype(float))
    inv_dp_scaled = minmax_scale(1.0 - diffpot_vals)  # invert: low diff pot = high risk
    stellate_scaled = minmax_scale(stellate_vals)
    velocity_scaled = minmax_scale(velocity_vals)

    cps = (CPS_W_S2 * s2_scaled
           + CPS_W_PF4 * pf4_scaled
           + CPS_W_DIFFPOT * inv_dp_scaled
           + CPS_W_STELLATE * stellate_scaled
           + CPS_W_VELOCITY * velocity_scaled)

    # CPS is NaN where any mandatory component is NaN
    # (S2 and P(F4) are the most restrictive at 1,254 samples)
    labels["cps"] = cps
    log.info("CPS: %d non-NaN values (mean=%.4f)", (~np.isnan(cps)).sum(), np.nanmean(cps))

    # Add metadata columns to labels
    labels["fibrosis_stage"] = meta["fibrosis_stage"].values
    labels["nas_score"] = meta["nas_score"].values
    labels["dataset"] = meta["dataset"].values
    labels["loco_fold_fibrosis"] = meta["loco_fold_fibrosis"].values

    # ------------------------------------------------------------------
    # Save outputs
    # ------------------------------------------------------------------
    log.info("Saving outputs to %s ...", OUT_DIR)

    # Feature matrix
    feat_path = os.path.join(OUT_DIR, "prognosis_feature_matrix.csv")
    feat_df.to_csv(feat_path)
    log.info("Saved: %s (%d x %d)", os.path.basename(feat_path), feat_df.shape[0], feat_df.shape[1])

    # Pseudo-labels
    labels_path = os.path.join(OUT_DIR, "prognosis_pseudo_labels.csv")
    labels.to_csv(labels_path)
    log.info("Saved: %s (%d x %d)", os.path.basename(labels_path), labels.shape[0], labels.shape[1])

    # Feature descriptions
    desc_df = pd.DataFrame(desc_rows, columns=["feature_name", "tier", "description"])
    desc_path = os.path.join(OUT_DIR, "prognosis_feature_description.csv")
    desc_df.to_csv(desc_path, index=False)
    log.info("Saved: %s (%d features)", os.path.basename(desc_path), len(desc_df))

    # HCC COLOC genes
    hcc_path = os.path.join(OUT_DIR, "hcc_coloc_genes.csv")
    if len(hcc_genes_out) > 0:
        hcc_genes_out.to_csv(hcc_path, index=False)
        log.info("Saved: %s (%d genes)", os.path.basename(hcc_path), len(hcc_genes_out))
    else:
        pd.DataFrame(columns=["gene", "PP.H4", "in_expression_matrix"]).to_csv(hcc_path, index=False)
        log.info("Saved: %s (empty — no HCC genes)", os.path.basename(hcc_path))

    elapsed = time.time() - t0
    log.info("Done in %.1f seconds. Feature matrix: %d samples x %d features, %d pseudo-labels",
             elapsed, n_samples, n_features, labels.shape[1])


if __name__ == "__main__":
    main()
