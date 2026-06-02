#!/usr/bin/env python3
"""
04b_scenic_grn_from_activity.py — Gene-activity-based GRN for hepatocytes (Module 3)

Builds a hepatocyte-specific gene regulatory network using SnapATAC2's gene
activity matrix as an RNA proxy. This approach is necessary because GSE244832
uses separate snATAC-seq (combinatorial indexing, 22bp barcodes) and snRNA-seq
(10x Chromium, 16bp barcodes) — barcodes cannot be matched cell-by-cell, so
native SCENIC+ (requiring paired ATAC+RNA per cell) cannot work.

Strategy:
  1. Load SnapATAC2 processed h5ad, subset to hepatocytes
  2. Compute gene activity matrix from AnnDataSet (fragment-based gene body
     accessibility) — a well-established approach (ArchR, Signac, SnapATAC2)
  3. Parse GENCODE v49 GTF for gene TSS positions
  4. Vectorized TF-enhancer correlations (rank matrix + matrix multiply)
  5. Enhancer-gene correlations for nearby gene-peak pairs (500kb window)
  6. Regulon assembly + differential activity (MASLD vs Normal)
  7. Export 4 CSVs matching L8 integration expectations

Inputs:
  - SnapATAC2 processed AnnData (88,814 cells with cell type + condition labels)
  - AnnDataSet (combined.h5ads) + per-donor backed h5ads for gene activity
  - GENCODE v49 GTF for TSS annotation

Outputs (in scenic_plus/):
  - hepatocyte_regulons.csv
  - enhancer_gene_links.csv
  - regulon_activity_scores.csv
  - disease_regulons.csv            (differential test now DONOR-level, +n_donors)
  - disease_regulons_excl_self.csv  (robustness: TF's own gene excluded; Track B)

NOTE (2026-05-30): The MASLD-vs-Normal differential regulon-activity test is run
at the DONOR level (18 donors; 5 NORMAL / 4 MASL + 9 MASH) via
utils_pseudobulk.aggregate_cells_to_donors + donor_groupwise_test. Per-cell
Mann-Whitney over thousands of cells from 18 donors was pseudoreplication
(review F033/F038). The per-cell activity matrix (regulon_activity_scores.csv) is
unchanged — only the aggregation+test changed. SCENIC regulon activity is also
dominated by the TF's own gene and targets were selected on the same cells being
tested (Track B circularity, F035/F039/F042); the self-excluded variant lets
reviewers gauge robustness.

Usage:
  python scripts/04b_scenic_grn_from_activity.py \
      --atac-h5ad results/snapatac2/snapatac2_processed.h5ad \
      --anndataset results/snapatac2/combined.h5ads \
      --gtf /gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz \
      --output-dir scenic_plus \
      --n-cells 15000

Author: MASLD-Atlas pipeline
"""

from __future__ import annotations

import argparse
import gc
import gzip
import logging
import os
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import mannwhitneyu as mwu, rankdata
from statsmodels.stats.multitest import multipletests

# Donor-level pseudobulk helpers (fixes per-cell pseudoreplication; F033/F038).
# Same dir as this script — make importable regardless of CWD.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils_pseudobulk import aggregate_cells_to_donors, donor_groupwise_test

log = logging.getLogger("scenic_grn")

# Curated donor metadata: donor_id (D01..D18, matches obs 'donor_id') -> condition
# (NORMAL / MASL / MASH). Used to build the donor->group map for the DONOR-level
# differential regulon-activity test (the donor is the experimental unit, n=18).
DONOR_META_PATH = (
    Path(__file__).resolve().parents[1] / "metadata" / "donor_metadata_curated.tsv"
)
DISEASE_CONDITIONS = ("MASL", "MASH")
HEALTHY_CONDITIONS = ("NORMAL",)


# ---------------------------------------------------------------------------
# Curated hepatocyte/liver TFs (66 TFs from Lambert et al. 2018 + liver TFs)
# ---------------------------------------------------------------------------
KNOWN_TFS = {
    "HNF4A", "HNF1A", "HNF1B", "HNF4G", "HNF6",
    "CEBPA", "CEBPB", "CEBPD", "CEBPG",
    "FOXA1", "FOXA2", "FOXA3",
    "PPARA", "PPARG", "PPARD",
    "RXRA", "RXRB", "RXRG",
    "NR1H4", "NR1H3", "NR1H2",  # FXR, LXR
    "SREBF1", "SREBF2",
    "MLX", "MLXIPL",
    "XBP1", "ATF4", "ATF6", "DDIT3",
    "STAT3", "STAT5A", "STAT5B",
    "JUN", "JUNB", "JUND", "FOS", "FOSB",
    "MYC", "MAX",
    "TP53", "TP63",
    "NFKB1", "NFKB2", "RELA", "RELB",
    "SP1", "SP3",
    "EGR1", "KLF4", "KLF6", "KLF15",
    "SOX9", "SOX4",
    "GATA4", "GATA6",
    "TCF7L2", "LEF1",
    "SMAD2", "SMAD3", "SMAD4",
    "ETS1", "ETS2", "ELF3",
    "IRF1", "IRF3",
    "ARNT", "AHR",
    "THRA", "THRB",  # Thyroid hormone receptors (resmetirom target)
    "ESR1", "AR",
    "RORA", "RORC",
}


# ---------------------------------------------------------------------------
# GTF parsing for TSS
# ---------------------------------------------------------------------------
def parse_tss_from_gtf(gtf_path: str) -> dict[str, tuple[str, int]]:
    """Parse gene-level TSS positions from GENCODE GTF.

    Returns dict: gene_name -> (chrom, tss_position)
    TSS = start for + strand genes, end for - strand genes.
    """
    log.info("Parsing TSS from GTF: %s", gtf_path)
    t0 = time.time()

    gene_tss = {}
    opener = gzip.open if gtf_path.endswith(".gz") else open

    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue

            chrom = fields[0]
            # Only standard chromosomes
            if not chrom.startswith("chr"):
                continue

            strand = fields[6]
            start = int(fields[3])  # 1-based
            end = int(fields[4])

            # Extract gene_name from attributes
            attrs = fields[8]
            gene_name = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_name"):
                    gene_name = attr.split('"')[1] if '"' in attr else attr.split(" ")[1]
                    break

            if gene_name is None:
                continue

            # TSS: start for + strand, end for - strand
            tss = start if strand == "+" else end
            gene_tss[gene_name] = (chrom, tss)

    log.info("  Parsed %d gene TSS positions in %.1fs", len(gene_tss), time.time() - t0)
    return gene_tss


# ---------------------------------------------------------------------------
# Peak coordinate parsing
# ---------------------------------------------------------------------------
def parse_peak_coords(peak_names: list[str]) -> pd.DataFrame:
    """Parse peak names like 'chr1:1000-2000' into a DataFrame."""
    records = []
    for i, pname in enumerate(peak_names):
        try:
            chrom, rest = pname.split(":")
            start, end = rest.split("-")
            records.append({"peak_idx": i, "chrom": chrom, "start": int(start), "end": int(end), "peak_name": pname})
        except Exception:
            pass
    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Build proximity index: peak -> nearby genes (within window)
# ---------------------------------------------------------------------------
def build_peak_gene_proximity(
    peak_df: pd.DataFrame,
    gene_tss: dict[str, tuple[str, int]],
    gene_names: list[str],
    window: int = 500_000,
) -> dict[int, list[tuple[str, int, int]]]:
    """Build index: peak_idx -> list of (gene_name, gene_col_idx, distance_to_tss).

    Only includes genes present in gene_names (from the activity matrix).
    """
    log.info("Building peak-gene proximity index (window=%dkb)...", window // 1000)
    t0 = time.time()

    # Build gene name -> column index mapping
    gene_name_set = set(gene_names)
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    # Group TSS by chromosome for efficient lookup
    chrom_genes = {}
    for gene, (chrom, tss) in gene_tss.items():
        if gene not in gene_name_set:
            continue
        chrom_genes.setdefault(chrom, []).append((gene, tss))

    # Sort by TSS for binary search
    for chrom in chrom_genes:
        chrom_genes[chrom].sort(key=lambda x: x[1])

    proximity = {}
    n_pairs = 0

    for _, row in peak_df.iterrows():
        chrom = row["chrom"]
        peak_mid = (row["start"] + row["end"]) // 2
        p_idx = row["peak_idx"]

        genes_on_chrom = chrom_genes.get(chrom, [])
        if not genes_on_chrom:
            continue

        nearby = []
        for gene, tss in genes_on_chrom:
            dist = abs(peak_mid - tss)
            if dist <= window:
                nearby.append((gene, gene_to_idx[gene], dist))

        if nearby:
            proximity[p_idx] = nearby
            n_pairs += len(nearby)

    log.info("  %d peak-gene pairs across %d peaks in %.1fs",
             n_pairs, len(proximity), time.time() - t0)
    return proximity


# ---------------------------------------------------------------------------
# Vectorized Spearman correlation
# ---------------------------------------------------------------------------
def _rank_matrix(X: np.ndarray) -> np.ndarray:
    """Rank each column of X (ties averaged). Returns float32 rank matrix."""
    n, m = X.shape
    R = np.empty((n, m), dtype=np.float32)
    for j in range(m):
        R[:, j] = rankdata(X[:, j], method="average")
    return R


def vectorized_spearman_cols(rank_A: np.ndarray, rank_B: np.ndarray) -> np.ndarray:
    """Spearman correlation between all columns of A and all columns of B.

    Returns (n_cols_A, n_cols_B) correlation matrix.
    rank_A and rank_B must already be rank-transformed.
    """
    n = rank_A.shape[0]
    # Center ranks
    A_c = rank_A - rank_A.mean(axis=0, keepdims=True)
    B_c = rank_B - rank_B.mean(axis=0, keepdims=True)

    # Norms
    A_norm = np.sqrt((A_c ** 2).sum(axis=0, keepdims=True))  # (1, n_cols_A)
    B_norm = np.sqrt((B_c ** 2).sum(axis=0, keepdims=True))  # (1, n_cols_B)

    # Avoid division by zero
    A_norm = np.maximum(A_norm, 1e-12)
    B_norm = np.maximum(B_norm, 1e-12)

    # Correlation matrix: (n_cols_A, n_cols_B)
    corr = (A_c / A_norm).T @ (B_c / B_norm)
    return corr


# ---------------------------------------------------------------------------
# Main GRN pipeline
# ---------------------------------------------------------------------------
def load_and_get_hepatocyte_indices(
    h5ad_path: str,
    n_cells: int = 15000,
    seed: int = 42,
) -> tuple[ad.AnnData, np.ndarray]:
    """Load SnapATAC2 processed h5ad, identify hepatocyte indices, subsample.

    Returns:
        (subsampled_adata, global_integer_indices)
        where global_integer_indices are positions in the FULL 88k-cell h5ad,
        used to subset the gene activity matrix by the same positions.
    """
    log.info("Loading ATAC data: %s", h5ad_path)
    adata = ad.read_h5ad(h5ad_path)
    log.info("  Full data: %d cells x %d peaks", adata.n_obs, adata.n_vars)

    # Get hepatocyte integer indices in the full h5ad
    hep_mask = adata.obs["cell_type"] == "Hepatocyte"
    hep_global_idx = np.where(hep_mask)[0]
    log.info("  Hepatocytes: %d cells", len(hep_global_idx))

    if len(hep_global_idx) > n_cells:
        # Stratified subsample across conditions
        rng = np.random.RandomState(seed)
        conditions = adata.obs["condition"].values[hep_global_idx]
        unique_conds = np.unique(conditions)

        # Allocate proportionally
        total = len(hep_global_idx)
        selected_local = []

        for c in unique_conds:
            c_local_idx = np.where(conditions == c)[0]
            n_sample = max(1, int(n_cells * len(c_local_idx) / total))
            n_sample = min(n_sample, len(c_local_idx))
            selected_local.extend(rng.choice(c_local_idx, n_sample, replace=False))

        selected_local = sorted(selected_local)[:n_cells]
        hep_global_idx = hep_global_idx[selected_local]
        log.info("  Subsampled to %d cells (stratified by condition)", len(hep_global_idx))

    adata_sub = adata[hep_global_idx].copy()
    log.info("  Condition distribution: %s",
             dict(adata_sub.obs["condition"].value_counts()))
    return adata_sub, hep_global_idx


def compute_gene_activity(
    anndataset_path: str,
    global_indices: np.ndarray,
) -> ad.AnnData:
    """Compute gene activity matrix from AnnDataSet using SnapATAC2.

    Uses positional integer indices (from the full ATAC h5ad) to subset cells,
    avoiding obs_name mismatch issues between the ATAC h5ad and AnnDataSet.

    Returns AnnData with gene activity scores (cells x genes).
    """
    import snapatac2 as snap

    log.info("Computing gene activity from AnnDataSet: %s", anndataset_path)
    t0 = time.time()

    dataset = snap.read_dataset(anndataset_path)
    log.info("  AnnDataSet: %d cells x %d features", dataset.n_obs, dataset.n_vars)

    gene_mat = snap.pp.make_gene_matrix(dataset, gene_anno=snap.genome.hg38)
    log.info("  Gene activity matrix: %d cells x %d genes", gene_mat.n_obs, gene_mat.n_vars)

    # Convert to in-memory if backed
    if hasattr(gene_mat, "to_memory"):
        gene_mat = gene_mat.to_memory()

    # Subset by positional integer indices (same row order as ATAC h5ad)
    log.info("  Subsetting gene activity to %d cells by position", len(global_indices))
    gene_mat_sub = gene_mat[global_indices].copy()

    # Assign unique obs_names to avoid downstream issues
    gene_mat_sub.obs_names = [f"cell_{i}" for i in range(gene_mat_sub.n_obs)]

    del gene_mat, dataset
    gc.collect()

    log.info("  Gene activity subset: %d cells x %d genes in %.1fs",
             gene_mat_sub.n_obs, gene_mat_sub.n_vars, time.time() - t0)
    return gene_mat_sub


def compute_tf_enhancer_correlations(
    atac_X_rank: np.ndarray,
    activity_X_rank: np.ndarray,
    tf_indices: dict[str, int],
    peak_names: list[str],
    top_n_peaks: int = 50000,
    peak_accessibility: np.ndarray | None = None,
    min_corr: float = 0.1,
) -> pd.DataFrame:
    """Vectorized TF-enhancer correlations.

    Correlates each TF's gene activity with each peak's accessibility
    using pre-ranked matrices and matrix multiplication.

    Returns DataFrame: tf_name, peak_name, peak_idx, tf_enhancer_corr
    """
    log.info("Computing TF-enhancer correlations (vectorized)...")
    t0 = time.time()

    n_cells = atac_X_rank.shape[0]
    n_peaks = atac_X_rank.shape[1]

    # Select top accessible peaks
    if peak_accessibility is not None and n_peaks > top_n_peaks:
        top_idx = np.argsort(peak_accessibility)[-top_n_peaks:]
        top_idx = np.sort(top_idx)
        atac_sub_rank = atac_X_rank[:, top_idx]
        log.info("  Using top %d peaks (of %d)", top_n_peaks, n_peaks)
    else:
        top_idx = np.arange(n_peaks)
        atac_sub_rank = atac_X_rank

    # Extract TF columns from activity rank matrix
    tf_names = list(tf_indices.keys())
    tf_col_idx = [tf_indices[tf] for tf in tf_names]
    tf_rank = activity_X_rank[:, tf_col_idx]  # (n_cells, n_tfs)

    # Filter out constant TFs
    tf_std = tf_rank.std(axis=0)
    valid_tf_mask = tf_std > 1e-10
    if valid_tf_mask.sum() < len(tf_names):
        log.info("  Removing %d constant TFs", (~valid_tf_mask).sum())
        tf_rank = tf_rank[:, valid_tf_mask]
        tf_names = [tf_names[i] for i in range(len(tf_names)) if valid_tf_mask[i]]

    # Filter out constant peaks
    peak_std = atac_sub_rank.std(axis=0)
    valid_peak_mask = peak_std > 1e-10
    atac_valid_rank = atac_sub_rank[:, valid_peak_mask]
    valid_peak_local_idx = np.where(valid_peak_mask)[0]
    log.info("  Valid: %d TFs x %d peaks", len(tf_names), valid_peak_local_idx.shape[0])

    # Vectorized Spearman: (n_tfs, n_valid_peaks)
    corr_mat = vectorized_spearman_cols(tf_rank, atac_valid_rank)

    # Extract significant correlations
    records = []
    for tf_i, tf_name in enumerate(tf_names):
        row = corr_mat[tf_i]
        sig_mask = np.abs(row) > min_corr
        sig_idx = np.where(sig_mask)[0]

        for s_idx in sig_idx:
            local_peak_idx = valid_peak_local_idx[s_idx]
            global_peak_idx = top_idx[local_peak_idx]
            records.append({
                "tf_name": tf_name,
                "peak_name": peak_names[global_peak_idx],
                "peak_idx": int(global_peak_idx),
                "tf_enhancer_corr": round(float(row[s_idx]), 4),
            })

    df = pd.DataFrame(records)
    log.info("  %d TF-enhancer links (|rho|>%.2f) in %.1fs",
             len(df), min_corr, time.time() - t0)
    return df


def compute_enhancer_gene_correlations(
    atac_X_rank: np.ndarray,
    activity_X_rank: np.ndarray,
    tf_enhancer_df: pd.DataFrame,
    proximity: dict[int, list[tuple[str, int, int]]],
    peak_df: pd.DataFrame,
    min_corr: float = 0.1,
    batch_size: int = 5000,
) -> pd.DataFrame:
    """Compute enhancer-gene correlations for peaks linked to TFs.

    For each TF-linked peak, correlate peak accessibility with nearby gene
    activity (proximity-based, within 500kb).

    Returns DataFrame with enhancer-gene link columns matching L8 expectations.
    """
    log.info("Computing enhancer-gene correlations...")
    t0 = time.time()

    # Get unique peaks from TF-enhancer results
    if tf_enhancer_df.empty:
        return pd.DataFrame()

    tf_peak_groups = tf_enhancer_df.groupby("peak_idx").agg({
        "tf_name": lambda x: list(x),
        "tf_enhancer_corr": lambda x: list(x),
    }).to_dict("index")

    # Peak coordinate lookup
    peak_coord_map = {}
    for _, row in peak_df.iterrows():
        peak_coord_map[row["peak_idx"]] = (row["chrom"], row["start"], row["end"])

    records = []
    n_computed = 0

    for p_idx, tf_info in tf_peak_groups.items():
        # Check if this peak has nearby genes
        if p_idx not in proximity:
            continue

        peak_acc_rank = atac_X_rank[:, p_idx]
        if peak_acc_rank.std() < 1e-10:
            continue

        chrom, start, end = peak_coord_map.get(p_idx, ("?", 0, 0))

        for gene_name, gene_col_idx, dist in proximity[p_idx]:
            gene_act_rank = activity_X_rank[:, gene_col_idx]
            if gene_act_rank.std() < 1e-10:
                continue

            # Spearman from pre-ranked data
            n = len(peak_acc_rank)
            d = peak_acc_rank - peak_acc_rank.mean()
            e = gene_act_rank - gene_act_rank.mean()
            denom = np.sqrt((d ** 2).sum() * (e ** 2).sum())
            rho = float((d * e).sum() / denom) if denom > 1e-12 else 0.0
            n_computed += 1

            if abs(rho) > min_corr:
                # Create one record per TF linked to this peak
                for tf_name, tf_corr in zip(tf_info["tf_name"], tf_info["tf_enhancer_corr"]):
                    records.append({
                        "enhancer_chr": chrom,
                        "enhancer_start": start,
                        "enhancer_end": end,
                        "target_gene": gene_name,
                        "tf_name": tf_name,
                        "correlation_rna_atac": round(rho, 4),
                        "tf_enhancer_corr": round(tf_corr, 4),
                        "distance_to_tss": dist,
                    })

    df = pd.DataFrame(records)
    log.info("  %d enhancer-gene links from %d correlations in %.1fs",
             len(df), n_computed, time.time() - t0)
    return df


def assemble_regulons(
    tf_enhancer_df: pd.DataFrame,
    enhancer_gene_df: pd.DataFrame,
) -> dict[str, dict]:
    """Assemble regulons: group peaks and target genes per TF."""
    log.info("Assembling regulons...")

    regulons = {}

    for tf_name in tf_enhancer_df["tf_name"].unique():
        tf_peaks = tf_enhancer_df.loc[
            tf_enhancer_df["tf_name"] == tf_name, "peak_name"
        ].unique().tolist()

        tf_targets = []
        if not enhancer_gene_df.empty:
            tf_targets = enhancer_gene_df.loc[
                enhancer_gene_df["tf_name"] == tf_name, "target_gene"
            ].unique().tolist()

        if tf_peaks:  # Keep TFs even without targets (they have enhancer evidence)
            regulons[tf_name] = {
                "peaks": tf_peaks,
                "target_genes": tf_targets,
            }

    log.info("  %d TFs with regulons, %d total target genes",
             len(regulons),
             sum(len(r["target_genes"]) for r in regulons.values()))
    return regulons


def build_donor_to_group(
    activity_adata: ad.AnnData,
    condition_col: str = "condition",
    donor_col: str = "donor_id",
) -> dict:
    """Map donor_id -> {'disease','healthy'} using the curated metadata.

    Authoritative source is metadata/donor_metadata_curated.tsv keyed on the
    'donor_id' column (D01..D18), which matches the activity obs 'donor_id'.
    NORMAL -> 'healthy'; MASL/MASH -> 'disease'. Falls back to the per-cell
    obs 'condition' (already row-aligned to the ATAC h5ad) for any donor not
    covered by the metadata, so the test never silently drops a donor.
    """
    donor_to_group: dict = {}

    # 1) Curated metadata (authoritative)
    if DONOR_META_PATH.exists():
        meta = pd.read_csv(DONOR_META_PATH, sep="\t", dtype=str)
        if {"donor_id", "condition"}.issubset(meta.columns):
            for _, r in meta.iterrows():
                cond = str(r["condition"]).strip().upper()
                if cond in DISEASE_CONDITIONS:
                    donor_to_group[str(r["donor_id"])] = "disease"
                elif cond in HEALTHY_CONDITIONS:
                    donor_to_group[str(r["donor_id"])] = "healthy"
            log.info("  donor->group from %s: %d donors mapped",
                     DONOR_META_PATH, len(donor_to_group))
        else:
            log.warning("  %s lacks donor_id/condition cols; using obs condition",
                        DONOR_META_PATH)
    else:
        log.warning("  Curated donor metadata not found at %s; "
                    "falling back to obs condition", DONOR_META_PATH)

    # 2) Fallback: derive from obs (donor_id, condition) for uncovered donors.
    if donor_col in activity_adata.obs and condition_col in activity_adata.obs:
        obs = activity_adata.obs[[donor_col, condition_col]].astype(str)
        for d, c in obs.drop_duplicates().itertuples(index=False):
            if d in donor_to_group:
                continue
            cu = c.strip().upper()
            if cu in DISEASE_CONDITIONS:
                donor_to_group[d] = "disease"
            elif cu in HEALTHY_CONDITIONS:
                donor_to_group[d] = "healthy"

    n_dis = sum(v == "disease" for v in donor_to_group.values())
    n_heal = sum(v == "healthy" for v in donor_to_group.values())
    log.info("  donor groups: %d disease, %d healthy", n_dis, n_heal)
    return donor_to_group


def _donor_differential(
    activity_matrix: np.ndarray,
    regulon_names: list[str],
    regulons: dict[str, dict],
    donor_array: np.ndarray,
    donor_to_group: dict,
    target_genes_per_regulon: dict[str, list[str]] | None = None,
) -> pd.DataFrame:
    """DONOR-level differential regulon activity (the donor is the unit).

    Aggregates the per-cell regulon-activity matrix to one value per donor
    (mean), then runs a per-regulon donor-level Mann-Whitney U (disease vs
    healthy) via the shared `donor_groupwise_test`. This replaces the previous
    per-CELL test, which was pseudoreplicated (~thousands of cells from 18
    donors). Returns a DataFrame with the SAME columns the downstream
    `disease_regulons.csv` schema expects, plus `n_donors`.

    Parameters
    ----------
    target_genes_per_regulon : optional override of each regulon's target-gene
        list used only for the reported `target_genes` / `n_target_genes`
        columns (used by the self-excluded robustness variant). The activity
        scores themselves are passed in via `activity_matrix`.
    """
    # cells x regulons -> donors x regulons (per-donor mean activity)
    donors, M_donor, n_cells_per_donor = aggregate_cells_to_donors(
        activity_matrix, donor_array, agg="mean", min_cells=1
    )
    log.info("  Aggregated %d cells -> %d donors for donor-level test",
             activity_matrix.shape[0], len(donors))

    dt = donor_groupwise_test(
        M_donor, donors, donor_to_group,
        group_a="disease", group_b="healthy",
        feature_names=regulon_names, min_per_group=2,
    )
    # donor_groupwise_test indexes rows positionally over feature_names, so the
    # row order matches regulon_names exactly.
    dt = dt.set_index("feature").reindex(regulon_names)
    n_donors_total = int(dt["n_donors_a"].iloc[0] + dt["n_donors_b"].iloc[0]) if len(dt) else 0

    records = []
    for tf_name in regulon_names:
        reg_data = regulons[tf_name]
        if target_genes_per_regulon is not None:
            target_genes = target_genes_per_regulon.get(tf_name, [])
        else:
            target_genes = reg_data.get("target_genes", [])
        n_enhancers = len(reg_data.get("peaks", []))
        row = dt.loc[tf_name]

        mean_masld = float(row["mean_disease"]) if pd.notna(row["mean_disease"]) else np.nan
        mean_normal = float(row["mean_healthy"]) if pd.notna(row["mean_healthy"]) else np.nan
        diff = float(row["mean_diff"]) if pd.notna(row["mean_diff"]) else np.nan
        pval = float(row["pvalue"]) if pd.notna(row["pvalue"]) else np.nan

        records.append({
            "regulon_id": f"{tf_name}_regulon",
            "tf_name": tf_name,
            "n_target_genes": len(target_genes),
            "target_genes": ";".join(target_genes[:200]),
            "n_enhancers": n_enhancers,
            "mean_activity_masld": round(mean_masld, 6) if not np.isnan(mean_masld) else np.nan,
            "mean_activity_normal": round(mean_normal, 6) if not np.isnan(mean_normal) else np.nan,
            "regulon_activity_diff": round(diff, 6) if not np.isnan(diff) else np.nan,
            "activity_pval": pval,
            # n_donors actually entering the donor-level test (disease + healthy)
            "n_donors": n_donors_total,
        })

    regulon_df = pd.DataFrame(records)

    # BH-FDR over the donor-level p-values (recomputed here so primary and
    # self-excluded variants each get an independent multiple-testing correction).
    if not regulon_df.empty and regulon_df["activity_pval"].notna().any():
        valid_mask = regulon_df["activity_pval"].notna()
        padj = np.full(len(regulon_df), np.nan)
        if valid_mask.sum() > 0:
            _, padj_valid, _, _ = multipletests(
                regulon_df.loc[valid_mask, "activity_pval"].values, method="fdr_bh"
            )
            padj[valid_mask.values] = padj_valid
        regulon_df["activity_padj"] = padj
    return regulon_df


def score_regulon_activity(
    activity_adata: ad.AnnData,
    regulons: dict[str, dict],
    condition_col: str = "condition",
    donor_col: str = "donor_id",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Score regulon activity per cell and test for differential activity.

    Uses AUCell-like mean expression of target genes as the per-cell activity
    score. The differential MASLD-vs-Normal test is run at the DONOR level
    (the donor, not the cell, is the experimental unit; 18 donors,
    5 NORMAL / 13 MASL+MASH) to avoid the pseudoreplication that previously
    declared nearly every regulon significant on thousands of cells
    (review F033/F038).

    Returns: (activity_df, disease_regulons_df, regulon_summary_df,
              disease_regulons_excl_self_df)
    """
    log.info("Scoring regulon activity per cell...")

    act_X = activity_adata.X
    if sp.issparse(act_X):
        act_X = act_X.toarray()

    gene_to_idx = {g: i for i, g in enumerate(activity_adata.var_names)}

    activity_matrix = np.zeros((activity_adata.shape[0], len(regulons)), dtype=np.float32)
    # Self-excluded activity matrix: same scores but the TF's OWN gene is dropped
    # from each regulon's target set before averaging. SCENIC regulon activity is
    # otherwise dominated by the TF's own promoter/gene, and targets were selected
    # on the same cells being tested (Track B circularity, review F035/F039/F042).
    activity_matrix_excl = np.zeros_like(activity_matrix)
    regulon_names: list[str] = []
    # target-gene list per regulon AFTER removing the TF's own gene (for reporting
    # in the self-excluded variant) + whether self-exclusion was actually applied.
    targets_excl_self: dict[str, list[str]] = {}
    tf_self_in_targets: dict[str, bool] = {}

    for reg_i, (tf_name, reg_data) in enumerate(regulons.items()):
        regulon_names.append(tf_name)
        target_genes = reg_data.get("target_genes", [])
        gene_indices = [gene_to_idx[g] for g in target_genes if g in gene_to_idx]

        # Drop the TF's own gene from the target set for the robustness variant.
        targets_no_self = [g for g in target_genes if g != tf_name]
        targets_excl_self[tf_name] = targets_no_self
        tf_self_in_targets[tf_name] = tf_name in target_genes
        gene_indices_excl = [gene_to_idx[g] for g in targets_no_self if g in gene_to_idx]

        if gene_indices:
            activity_matrix[:, reg_i] = act_X[:, gene_indices].mean(axis=1)
        if gene_indices_excl:
            activity_matrix_excl[:, reg_i] = act_X[:, gene_indices_excl].mean(axis=1)

    activity_df = pd.DataFrame(
        activity_matrix,
        index=activity_adata.obs_names,
        columns=regulon_names,
    )

    # ------------------------------------------------------------------
    # DONOR-level differential regulon activity (MASLD vs Normal).
    # The unit is the donor (18), NOT the cell. Aggregate per-cell scores to
    # per-donor means, then Mann-Whitney across donor groups.
    # ------------------------------------------------------------------
    log.info("Testing differential regulon activity at the DONOR level...")
    donor_to_group = build_donor_to_group(activity_adata, condition_col, donor_col)
    if donor_col in activity_adata.obs:
        donor_array = activity_adata.obs[donor_col].astype(str).to_numpy()
    else:
        # No donor labels available — degrade gracefully to a single pseudo-donor
        # per condition so the test still runs (but warn loudly: this is NOT a
        # valid donor-level test).
        log.warning("  obs lacks '%s'; falling back to condition as pseudo-donor "
                    "(NOT a valid donor-level test)", donor_col)
        donor_array = activity_adata.obs.get(
            condition_col, pd.Series("Unknown", index=activity_adata.obs_names)
        ).astype(str).to_numpy()

    regulon_df = _donor_differential(
        activity_matrix, regulon_names, regulons, donor_array, donor_to_group,
    )

    # Self-excluded robustness variant (secondary output; primary is unchanged).
    regulon_df_excl = _donor_differential(
        activity_matrix_excl, regulon_names, regulons, donor_array, donor_to_group,
        target_genes_per_regulon=targets_excl_self,
    )
    # Flag regulons whose self-exclusion flips significance or sign, so reviewers
    # can see which calls are driven by the TF's own gene (TF-self-dominated).
    if not regulon_df.empty:
        primary_sig = regulon_df.set_index("tf_name")["activity_padj"] < 0.05
        excl_sig = regulon_df_excl.set_index("tf_name")["activity_padj"] < 0.05
        primary_dir = np.sign(regulon_df.set_index("tf_name")["regulon_activity_diff"])
        excl_dir = np.sign(regulon_df_excl.set_index("tf_name")["regulon_activity_diff"])
        dominated = {}
        for tf in regulon_df["tf_name"]:
            has_self = tf_self_in_targets.get(tf, False)
            sig_flip = bool(primary_sig.get(tf, False)) and not bool(excl_sig.get(tf, False))
            dir_flip = pd.notna(primary_dir.get(tf)) and pd.notna(excl_dir.get(tf)) \
                and primary_dir.get(tf) != excl_dir.get(tf)
            dominated[tf] = bool(has_self and (sig_flip or dir_flip))
        regulon_df["is_tf_self_dominated"] = regulon_df["tf_name"].map(dominated)

    # Expand regulon_df: one row per (tf, target_gene) for L8 consumption
    expanded_rows = []
    for _, row in regulon_df.iterrows():
        tf_name = row["tf_name"]
        targets = row["target_genes"].split(";") if row["target_genes"] else []
        for target in targets:
            target = target.strip()
            if target:
                expanded_rows.append({
                    "tf_name": tf_name,
                    "target_gene": target,
                    "regulon_id": row["regulon_id"],
                    "n_target_genes": row["n_target_genes"],
                    "target_genes": row["target_genes"],
                    "n_enhancers": row["n_enhancers"],
                    "mean_activity_masld": row["mean_activity_masld"],
                    "mean_activity_normal": row["mean_activity_normal"],
                    "regulon_activity_diff": row["regulon_activity_diff"],
                    "activity_pval": row["activity_pval"],
                    "activity_padj": row.get("activity_padj", np.nan),
                    "n_donors": row.get("n_donors", np.nan),
                })

    expanded_df = pd.DataFrame(expanded_rows)

    # Disease regulons (significant at the DONOR level)
    def _significant(df: pd.DataFrame) -> pd.DataFrame:
        if df.empty or "activity_padj" not in df.columns:
            return df.iloc[0:0].copy()
        out = df[df["activity_padj"] < 0.05].copy()
        return out.sort_values("activity_padj").reset_index(drop=True)

    disease_df = _significant(regulon_df)
    disease_df_excl = _significant(regulon_df_excl)

    log.info("Regulon activity: %d regulons scored, %d disease-associated "
             "(donor-level padj<0.05); self-excluded: %d",
             len(regulon_names), len(disease_df), len(disease_df_excl))

    return activity_df, disease_df, expanded_df, disease_df_excl


def export_results(
    out_dir: str,
    regulon_df: pd.DataFrame,
    enhancer_df: pd.DataFrame,
    activity_df: pd.DataFrame,
    disease_df: pd.DataFrame,
    disease_df_excl: pd.DataFrame | None = None,
):
    """Save the canonical L8 output CSVs (paths/columns preserved).

    Primary outputs are unchanged in path and schema; `disease_regulons.csv`
    now additionally carries `n_donors`. An OPTIONAL secondary file
    `disease_regulons_excl_self.csv` (TF's own gene removed from each regulon)
    is written for the Track-B circularity robustness check.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # 1. Hepatocyte regulons (one row per TF-target pair)
    path = out / "hepatocyte_regulons.csv"
    regulon_df.to_csv(path, index=False)
    log.info("  Saved: %s (%d rows)", path, len(regulon_df))

    # 2. Enhancer-gene links
    path = out / "enhancer_gene_links.csv"
    if enhancer_df is not None and not enhancer_df.empty:
        enhancer_df.to_csv(path, index=False)
    else:
        pd.DataFrame(columns=[
            "enhancer_chr", "enhancer_start", "enhancer_end",
            "target_gene", "tf_name", "correlation_rna_atac",
            "tf_enhancer_corr", "distance_to_tss",
        ]).to_csv(path, index=False)
    log.info("  Saved: %s (%d links)", path,
             len(enhancer_df) if enhancer_df is not None else 0)

    # 3. Regulon activity scores (cells x regulons)
    path = out / "regulon_activity_scores.csv"
    activity_df.to_csv(path)
    log.info("  Saved: %s (%s)", path, activity_df.shape)

    # 4. Disease regulons (donor-level padj < 0.05 subset; now incl. n_donors)
    path = out / "disease_regulons.csv"
    disease_df.to_csv(path, index=False)
    log.info("  Saved: %s (%d disease regulons)", path, len(disease_df))

    # 5. (optional) Self-excluded robustness variant — TF's own gene removed from
    #    each regulon before the donor-level test. Lets reviewers gauge how much
    #    each disease call depends on the TF's own promoter/gene (Track B).
    if disease_df_excl is not None:
        path = out / "disease_regulons_excl_self.csv"
        disease_df_excl.to_csv(path, index=False)
        log.info("  Saved: %s (%d disease regulons, TF-self excluded)",
                 path, len(disease_df_excl))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Gene-activity-based GRN for hepatocytes (Module 3)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--atac-h5ad", required=True,
        help="SnapATAC2 processed h5ad with cell type + condition labels",
    )
    parser.add_argument(
        "--anndataset", required=True,
        help="AnnDataSet path (combined.h5ads) for gene activity computation",
    )
    parser.add_argument(
        "--gtf", required=True,
        help="GENCODE GTF (gzipped) for TSS annotation",
    )
    parser.add_argument(
        "--output-dir", default="scenic_plus",
        help="Output directory for 4 CSV files",
    )
    parser.add_argument(
        "--n-cells", type=int, default=15000,
        help="Max hepatocytes to subsample (balanced across conditions)",
    )
    parser.add_argument(
        "--window", type=int, default=500000,
        help="Peak-gene proximity window in bp",
    )
    parser.add_argument(
        "--min-corr", type=float, default=0.1,
        help="Min |Spearman rho| for significant correlations",
    )
    parser.add_argument(
        "--top-peaks", type=int, default=50000,
        help="Top accessible peaks to use for TF-enhancer correlations",
    )
    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stdout),
        ],
    )

    log.info("=" * 70)
    log.info("Module 3: Gene-Activity-Based Hepatocyte GRN")
    log.info("=" * 70)
    t_start = time.time()

    # -----------------------------------------------------------------------
    # Step 1: Load ATAC data, subset to hepatocytes (returns global indices)
    # -----------------------------------------------------------------------
    atac_adata, hep_global_idx = load_and_get_hepatocyte_indices(
        args.atac_h5ad, n_cells=args.n_cells,
    )

    # -----------------------------------------------------------------------
    # Step 2: Compute gene activity matrix from AnnDataSet
    #   Uses same positional indices to avoid obs_name mismatch
    # -----------------------------------------------------------------------
    activity_adata = compute_gene_activity(args.anndataset, hep_global_idx)

    # Align obs_names and transfer condition + donor labels (row-aligned by the
    # positional subsetting above, hence the n_obs assert). donor_id is required
    # for the DONOR-level differential test; without it the test degrades to a
    # pseudoreplicated per-condition comparison.
    assert atac_adata.n_obs == activity_adata.n_obs, (
        f"Cell count mismatch: ATAC={atac_adata.n_obs}, activity={activity_adata.n_obs}"
    )
    activity_adata.obs_names = atac_adata.obs_names.copy()
    activity_adata.obs["condition"] = atac_adata.obs["condition"].values
    if "donor_id" in atac_adata.obs.columns:
        activity_adata.obs["donor_id"] = atac_adata.obs["donor_id"].values
        log.info("  Transferred donor_id (%d donors) to activity adata",
                 atac_adata.obs["donor_id"].astype(str).nunique())
    else:
        log.warning("  ATAC obs has no 'donor_id'; donor-level test will degrade "
                    "to per-condition (pseudoreplicated). Check the input h5ad.")
    log.info("Aligned %d cells between ATAC and gene activity", atac_adata.n_obs)

    # -----------------------------------------------------------------------
    # Step 3: Parse TSS from GTF
    # -----------------------------------------------------------------------
    gene_tss = parse_tss_from_gtf(args.gtf)

    # -----------------------------------------------------------------------
    # Step 4: Identify TFs present in gene activity matrix
    # -----------------------------------------------------------------------
    activity_genes = set(activity_adata.var_names)
    tfs_present = sorted(KNOWN_TFS & activity_genes)
    tf_to_idx = {tf: list(activity_adata.var_names).index(tf) for tf in tfs_present}
    log.info("TFs present in gene activity: %d / %d", len(tfs_present), len(KNOWN_TFS))
    log.info("  TFs: %s", ", ".join(tfs_present))

    if len(tfs_present) == 0:
        log.error("No known TFs found in gene activity matrix. Exiting.")
        sys.exit(1)

    # -----------------------------------------------------------------------
    # Step 5: Select top accessible peaks from sparse matrix, then densify
    #   Full ATAC matrix is ~15K x 6M peaks = 339 GiB dense; select top N first
    # -----------------------------------------------------------------------
    log.info("Selecting top %d accessible peaks from sparse ATAC matrix...", args.top_peaks)
    t_rank = time.time()

    atac_sparse = atac_adata.X
    if not sp.issparse(atac_sparse):
        atac_sparse = sp.csr_matrix(atac_sparse)
    else:
        atac_sparse = sp.csr_matrix(atac_sparse)  # ensure CSR for efficient column ops

    # Compute mean accessibility per peak from sparse matrix
    peak_accessibility = np.asarray(atac_sparse.mean(axis=0)).ravel()
    top_peak_idx = np.argsort(peak_accessibility)[-args.top_peaks:]
    top_peak_idx = np.sort(top_peak_idx)

    all_peak_names = list(atac_adata.var_names)
    top_peak_names = [all_peak_names[i] for i in top_peak_idx]

    log.info("  Selected %d peaks (mean accessibility range: %.4f - %.4f)",
             len(top_peak_idx),
             peak_accessibility[top_peak_idx[0]],
             peak_accessibility[top_peak_idx[-1]])

    # Build global-to-local index for top peaks
    global_to_local_peak = {g: l for l, g in enumerate(top_peak_idx)}

    # Densify only the top peaks subset
    atac_X = np.asarray(atac_sparse[:, top_peak_idx].toarray(), dtype=np.float32)
    log.info("  Dense ATAC subset: %d x %d (%.1f GiB)",
             *atac_X.shape, atac_X.nbytes / 1e9)

    del atac_sparse
    gc.collect()

    # Densify gene activity matrix
    act_X = activity_adata.X
    if sp.issparse(act_X):
        act_X = act_X.toarray()
    act_X = np.asarray(act_X, dtype=np.float32)

    # Rank ATAC peaks (top subset only)
    log.info("  Ranking ATAC matrix (%d x %d)...", *atac_X.shape)
    atac_X_rank = _rank_matrix(atac_X)

    # Rank gene activity
    log.info("  Ranking activity matrix (%d x %d)...", *act_X.shape)
    activity_X_rank = _rank_matrix(act_X)

    del atac_X, act_X
    gc.collect()
    log.info("  Rank matrices computed in %.1fs", time.time() - t_rank)

    # -----------------------------------------------------------------------
    # Step 6: Parse peak coordinates + build proximity index
    #   Only for the top peaks subset
    # -----------------------------------------------------------------------
    peak_df = parse_peak_coords(top_peak_names)
    proximity = build_peak_gene_proximity(
        peak_df, gene_tss,
        gene_names=list(activity_adata.var_names),
        window=args.window,
    )

    # -----------------------------------------------------------------------
    # Step 7: Vectorized TF-enhancer correlations
    #   Now operates on the pre-selected top peaks (no further subsetting)
    # -----------------------------------------------------------------------
    tf_enhancer_df = compute_tf_enhancer_correlations(
        atac_X_rank=atac_X_rank,
        activity_X_rank=activity_X_rank,
        tf_indices=tf_to_idx,
        peak_names=top_peak_names,
        top_n_peaks=len(top_peak_names),  # Already pre-selected
        peak_accessibility=None,           # No further filtering needed
        min_corr=args.min_corr,
    )

    # -----------------------------------------------------------------------
    # Step 8: Enhancer-gene correlations (proximity-constrained)
    # -----------------------------------------------------------------------
    enhancer_gene_df = compute_enhancer_gene_correlations(
        atac_X_rank=atac_X_rank,
        activity_X_rank=activity_X_rank,
        tf_enhancer_df=tf_enhancer_df,
        proximity=proximity,
        peak_df=peak_df,
        min_corr=args.min_corr,
    )

    # -----------------------------------------------------------------------
    # Step 9: Assemble regulons
    # -----------------------------------------------------------------------
    regulons = assemble_regulons(tf_enhancer_df, enhancer_gene_df)

    if len(regulons) < 10:
        log.warning("Only %d TFs with regulons (threshold: 10). Consider lowering --min-corr.", len(regulons))

    # -----------------------------------------------------------------------
    # Step 10: Differential regulon activity (DONOR-level; + self-excluded arm)
    # -----------------------------------------------------------------------
    activity_df, disease_df, regulon_summary_df, disease_df_excl = score_regulon_activity(
        activity_adata, regulons,
    )

    # -----------------------------------------------------------------------
    # Step 11: Export results
    # -----------------------------------------------------------------------
    export_results(
        args.output_dir,
        regulon_df=regulon_summary_df,
        enhancer_df=enhancer_gene_df,
        activity_df=activity_df,
        disease_df=disease_df,
        disease_df_excl=disease_df_excl,
    )

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    elapsed = time.time() - t_start
    log.info("=" * 70)
    log.info("Module 3 GRN complete in %.1f minutes", elapsed / 60)
    log.info("  TFs analyzed: %d", len(tfs_present))
    log.info("  TFs with regulons: %d", len(regulons))
    log.info("  TF-enhancer links: %d", len(tf_enhancer_df))
    log.info("  Enhancer-gene links: %d", len(enhancer_gene_df))
    log.info("  Unique target genes: %d",
             enhancer_gene_df["target_gene"].nunique() if not enhancer_gene_df.empty else 0)
    log.info("  Disease regulons (padj<0.05): %d", len(disease_df))
    log.info("  Output directory: %s", args.output_dir)
    log.info("=" * 70)

    # Biological sanity checks
    expected_tfs = {"HNF4A", "FOXA1", "FOXA2", "FOXA3", "CEBPA", "CEBPB", "THRB", "PPARA"}
    found_expected = expected_tfs & set(regulons.keys())
    missing_expected = expected_tfs - set(regulons.keys())
    log.info("Biological validation:")
    log.info("  Expected hepatocyte TFs with regulons: %s", ", ".join(sorted(found_expected)))
    if missing_expected:
        log.warning("  Missing expected TFs: %s", ", ".join(sorted(missing_expected)))


if __name__ == "__main__":
    main()
