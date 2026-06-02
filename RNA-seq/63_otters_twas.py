#!/usr/bin/env python3
"""
63_otters_twas.py
Phase 2A-Step3: Run TWAS using trained OTTERS weights against GWAS summary statistics.

For each gene, computes a TWAS z-score:
  z_twas = w' * z_gwas / sqrt(w' * R * w)
where w = prediction weights, z_gwas = GWAS Z-scores at weight SNPs, R = LD matrix.

Methods combined via ACAT (Cauchy combination test) across P+T and lassosum.

Input:
  - data/broadaway_eqtl/otters_weights/chr{1-22}/*.weights.tsv  (from Script 62)
  - data/broadaway_eqtl/otters_format/gene_anno.txt
  - GWAS/MR_Data/hg19/*.tsv.gz                                  (GWAS hg19 liftover)
  - data/1kg_eur/chr{1-22}_eur.{pgen,psam,pvar}                 (1KG EUR LD)

Output:
  - RNA-seq/results/causal_inference/otters_broadaway/{gwas}/otters_twas_combined.csv

Usage:
  GWAS_NAME=ghodsian python 63_otters_twas.py
  GWAS_NAME=ukbb_alt python 63_otters_twas.py
  # Or via SLURM:
  GWAS_NAME=ghodsian sbatch run_otters_twas.sbatch

Requires: numpy, pandas, scipy, subprocess (plink2 via module load)

Author: James Lee
Date: 2026-03-24
"""

import os
import sys
import time
import json
import argparse
import subprocess
import tempfile
import shutil
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import norm, cauchy


# =============================================================================
# Configuration
# =============================================================================

MIN_SNPS_FOR_TWAS = 3  # minimum overlapping SNPs to compute TWAS
MIN_WEIGHT_NORM = 1e-10  # minimum weight vector norm (avoid near-zero denominators)
CHECKPOINT_INTERVAL = 500
PROGRESS_INTERVAL = 100

# GWAS file mapping: name -> (filename, column mapping, sample size)
# Column mapping: {standard_name: gwas_column_name}
# Standard columns needed: chr, pos, effect_allele, other_allele, beta, se, pval
GWAS_REGISTRY = {
    "ghodsian": {
        "file": "Ghodsian_2021_NAFLD_harmonised_hg19.tsv.gz",
        "cols": {
            "chr": "hg19_chr",
            "pos": "hg19_pos",
            "effect_allele": "hm_effect_allele",
            "other_allele": "hm_other_allele",
            "beta": "hm_beta",
            "se": "standard_error",
            "pval": "p_value",
        },
        "N": 778614,
        "type": "cc",  # case-control
    },
    "ukbb_alt": {
        "file": "GCST90019492_UKBB_ALT_harmonised_hg19.tsv.gz",
        "cols": {
            "chr": "hg19_chr",
            "pos": "hg19_pos",
            "effect_allele": "effect_allele",
            "other_allele": "other_allele",
            "beta": "beta",
            "se": "standard_error",
            "pval": "p_value",
        },
        "N": 343850,
        "type": "quant",
    },
    "ukbb_ast": {
        "file": "GCST90019497_UKBB_AST_harmonised_hg19.tsv.gz",
        "cols": {
            "chr": "hg19_chr",
            "pos": "hg19_pos",
            "effect_allele": "effect_allele",
            "other_allele": "other_allele",
            "beta": "beta",
            "se": "standard_error",
            "pval": "p_value",
        },
        "N": 343850,
        "type": "quant",
    },
    "ukbb_ggt": {
        "file": "GCST90019507_UKBB_GGT_harmonised_hg19.tsv.gz",
        "cols": {
            "chr": "hg19_chr",
            "pos": "hg19_pos",
            "effect_allele": "effect_allele",
            "other_allele": "other_allele",
            "beta": "beta",
            "se": "standard_error",
            "pval": "p_value",
        },
        "N": 343850,
        "type": "quant",
    },
}


def get_project_root():
    """Get project root from env or default."""
    return os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
    )


# =============================================================================
# GWAS loading
# =============================================================================

def load_gwas(gwas_path, col_map, chrom=None):
    """Load GWAS summary statistics with hg19 coordinates.

    Args:
        gwas_path: Path to GWAS file (gzipped TSV)
        col_map: Dictionary mapping standard names to GWAS column names
        chrom: Optional chromosome filter (int)

    Returns:
        DataFrame with columns: chr, pos, effect_allele, other_allele, beta, se, pval, zscore
    """
    # Determine which columns to load
    usecols = list(col_map.values())
    df = pd.read_csv(gwas_path, sep="\t", usecols=usecols)

    # Rename to standard names
    inv_map = {v: k for k, v in col_map.items()}
    df = df.rename(columns=inv_map)

    # Coerce types
    df["chr"] = pd.to_numeric(df["chr"], errors="coerce")
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df["beta"] = pd.to_numeric(df["beta"], errors="coerce")
    df["se"] = pd.to_numeric(df["se"], errors="coerce")
    df["pval"] = pd.to_numeric(df["pval"], errors="coerce")

    # Drop rows with missing critical values
    df = df.dropna(subset=["chr", "pos", "beta", "se"])

    # Filter to valid entries
    df = df[df["se"] > 0].copy()

    # Compute Z-score
    df["zscore"] = df["beta"] / df["se"]

    # Uppercase alleles
    df["effect_allele"] = df["effect_allele"].astype(str).str.upper()
    df["other_allele"] = df["other_allele"].astype(str).str.upper()

    # Chromosome filter
    if chrom is not None:
        df = df[df["chr"] == chrom].copy()

    # Convert pos to int
    df["pos"] = df["pos"].astype(np.int64)
    df["chr"] = df["chr"].astype(int)

    return df


# =============================================================================
# LD computation (same as Script 62)
# =============================================================================

def compute_ld_matrix(pgen_prefix, snp_ids, tmpdir):
    """Compute LD r matrix for a set of SNPs using plink2.

    Returns:
        ld_matrix: np.ndarray of shape (n_snps, n_snps) — Pearson r
        kept_ids: List of SNP IDs retained
    """
    if len(snp_ids) < 2:
        return np.array([[1.0]]), list(snp_ids)

    # Write SNP list
    snp_file = os.path.join(tmpdir, "snps.txt")
    with open(snp_file, "w") as f:
        for sid in snp_ids:
            f.write(sid + "\n")

    out_prefix = os.path.join(tmpdir, "ld_out")

    cmd = (
        f"module load PLINK/2.00a3.7-gfbf-2023a && "
        f"plink2 --pfile {pgen_prefix} "
        f"--extract {snp_file} "
        f"--r-unphased square "
        f"--out {out_prefix} "
        f"--threads 1 --memory 4000"
    )
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, executable="/bin/bash"
    )

    vcor_file = None
    for suffix in [".unphased.vcor", ".unphased.vcor1"]:
        candidate = out_prefix + suffix
        if os.path.exists(candidate):
            vcor_file = candidate
            break

    if vcor_file is None:
        # Fallback for older plink2
        cmd_fallback = (
            f"module load PLINK/2.00a3.7-gfbf-2023a && "
            f"plink2 --pfile {pgen_prefix} "
            f"--extract {snp_file} "
            f"--r square "
            f"--out {out_prefix} "
            f"--threads 1 --memory 4000"
        )
        result = subprocess.run(
            cmd_fallback, shell=True, capture_output=True, text=True,
            executable="/bin/bash"
        )
        for suffix in [".vcor", ".vcor1", ".ld"]:
            candidate = out_prefix + suffix
            if os.path.exists(candidate):
                vcor_file = candidate
                break

    if vcor_file is None:
        raise RuntimeError(
            f"plink2 LD failed. stderr:\n{result.stderr[:500]}"
        )

    ld_matrix = np.loadtxt(vcor_file, delimiter="\t")

    vars_file = None
    for suffix in [".unphased.vcor.vars", ".unphased.vcor1.vars",
                   ".vcor.vars", ".vcor1.vars", ".ld.vars"]:
        candidate = out_prefix + suffix
        if os.path.exists(candidate):
            vars_file = candidate
            break

    if vars_file and os.path.exists(vars_file):
        with open(vars_file) as f:
            kept_ids = [line.strip() for line in f if line.strip()]
    else:
        kept_ids = list(snp_ids)
        if ld_matrix.shape[0] != len(kept_ids):
            raise RuntimeError(
                f"LD matrix dim {ld_matrix.shape} != SNP count {len(kept_ids)}"
            )

    return ld_matrix, kept_ids


# =============================================================================
# TWAS z-score computation
# =============================================================================

def compute_twas_zscore(weights, gwas_z, ld_r):
    """Compute TWAS z-score.

    z_twas = w' * z_gwas / sqrt(w' * R * w)

    Args:
        weights: np.ndarray of prediction weights (n_snps,)
        gwas_z: np.ndarray of GWAS z-scores (n_snps,)
        ld_r: np.ndarray LD correlation matrix (n_snps, n_snps)

    Returns:
        z_twas: float (TWAS z-score) or None if denominator is too small
    """
    numerator = np.dot(weights, gwas_z)
    denominator_sq = np.dot(weights, ld_r @ weights)

    if denominator_sq <= MIN_WEIGHT_NORM:
        return None

    denominator = np.sqrt(denominator_sq)
    z_twas = numerator / denominator

    return z_twas


def acat_combine(pvals):
    """ACAT (Aggregated Cauchy Association Test) for combining p-values.

    Liu & Xie (2020) JASA. Robust to correlated tests.

    T = sum(tan((0.5 - p_i) * pi)) / k
    p_acat = 0.5 - arctan(T) / pi

    Args:
        pvals: list/array of p-values (removes NaN and 0/1 extremes)

    Returns:
        combined_pval: float
    """
    pvals = np.array([p for p in pvals if p is not None and np.isfinite(p)])
    if len(pvals) == 0:
        return None

    # Clamp extreme values to avoid tan(±inf)
    pvals = np.clip(pvals, 1e-300, 1 - 1e-15)

    if len(pvals) == 1:
        return float(pvals[0])

    # ACAT statistic
    T = np.mean(np.tan((0.5 - pvals) * np.pi))

    # Combined p-value via Cauchy CDF
    p_acat = 0.5 - np.arctan(T) / np.pi

    return float(np.clip(p_acat, 1e-300, 1.0))


# =============================================================================
# Gene-level TWAS
# =============================================================================

def load_pvar_positions(pvar_path):
    """Load pvar file for position/allele lookups."""
    df = pd.read_csv(
        pvar_path, sep="\t", comment="#", header=None,
        names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"],
        usecols=[1, 2, 3, 4],
    )
    df.columns = ["POS", "ID", "REF", "ALT"]
    df["POS"] = df["POS"].astype(np.int64)
    snp_mask = (df["REF"].str.len() == 1) & (df["ALT"].str.len() == 1)
    df = df[snp_mask].copy()
    return df


def run_twas_for_gene(weight_path, gwas_chr_df, pvar_df, pgen_prefix, tmpdir):
    """Run TWAS for a single gene.

    Args:
        weight_path: Path to weight file (SNP_ID, POS, A1, A2, Weight_PT, Weight_lassosum)
        gwas_chr_df: GWAS DataFrame for this chromosome
        pvar_df: pvar DataFrame for position lookups
        pgen_prefix: pgen file prefix for LD computation
        tmpdir: temp directory

    Returns:
        dict with TWAS results or None if skipped
    """
    # Load weights
    weights_df = pd.read_csv(weight_path, sep="\t")
    if len(weights_df) < MIN_SNPS_FOR_TWAS:
        return None

    # Merge weights with GWAS on position
    # Weight SNP_ID is in 1KG format (e.g., '21:9411618:C:T')
    # Weight A1 = effect allele (1KG ALT), A2 = other allele (1KG REF)
    # GWAS has pos, effect_allele, other_allele, zscore

    merged = weights_df.merge(gwas_chr_df, left_on="POS", right_on="pos", how="inner")

    if len(merged) < MIN_SNPS_FOR_TWAS:
        return None

    # Allele alignment: ensure GWAS effect allele matches weight A1
    # Case 1: GWAS effect_allele == Weight A1 and GWAS other_allele == Weight A2 (concordant)
    concordant = (
        (merged["effect_allele"] == merged["A1"]) &
        (merged["other_allele"] == merged["A2"])
    )
    # Case 2: GWAS effect_allele == Weight A2 and GWAS other_allele == Weight A1 (flipped)
    flipped = (
        (merged["effect_allele"] == merged["A2"]) &
        (merged["other_allele"] == merged["A1"])
    )

    aligned = merged[concordant | flipped].copy()
    if len(aligned) < MIN_SNPS_FOR_TWAS:
        return None

    # Flip GWAS Z-score where alleles are reversed
    flip_mask = flipped[concordant | flipped].values
    aligned.loc[flip_mask, "zscore"] = -aligned.loc[flip_mask, "zscore"]

    # Deduplicate by position
    aligned = aligned.drop_duplicates(subset=["POS"], keep="first")

    if len(aligned) < MIN_SNPS_FOR_TWAS:
        return None

    n_snps = len(aligned)
    snp_ids = aligned["SNP_ID"].tolist()
    gwas_z = aligned["zscore"].values.astype(np.float64)
    w_pt = aligned["Weight_PT"].values.astype(np.float64)
    w_lasso = aligned["Weight_lassosum"].values.astype(np.float64)

    # Compute LD matrix
    # Clean tmpdir
    for f in os.listdir(tmpdir):
        try:
            os.remove(os.path.join(tmpdir, f))
        except OSError:
            pass

    try:
        ld_r, kept_ids = compute_ld_matrix(pgen_prefix, snp_ids, tmpdir)
    except RuntimeError:
        return None

    # Subset to kept SNPs
    if len(kept_ids) != len(snp_ids):
        kept_set = set(kept_ids)
        keep_mask = aligned["SNP_ID"].isin(kept_set)
        aligned = aligned[keep_mask].copy()
        gwas_z = aligned["zscore"].values.astype(np.float64)
        w_pt = aligned["Weight_PT"].values.astype(np.float64)
        w_lasso = aligned["Weight_lassosum"].values.astype(np.float64)
        n_snps = len(aligned)

    if n_snps < MIN_SNPS_FOR_TWAS:
        return None

    # Sanity check dimensions
    if ld_r.shape[0] != n_snps:
        return None

    # Fix LD matrix
    if np.any(~np.isfinite(ld_r)):
        ld_r = np.nan_to_num(ld_r, nan=0.0, posinf=1.0, neginf=-1.0)
    np.fill_diagonal(ld_r, 1.0)

    # Compute TWAS z-scores for each method
    results = {"n_snps": n_snps}

    # P+T
    if np.any(w_pt != 0):
        z_pt = compute_twas_zscore(w_pt, gwas_z, ld_r)
        if z_pt is not None and np.isfinite(z_pt):
            results["pt_z"] = z_pt
            results["pt_pval"] = 2 * norm.sf(abs(z_pt))
        else:
            results["pt_z"] = np.nan
            results["pt_pval"] = np.nan
    else:
        results["pt_z"] = np.nan
        results["pt_pval"] = np.nan

    # lassosum
    if np.any(w_lasso != 0):
        z_lasso = compute_twas_zscore(w_lasso, gwas_z, ld_r)
        if z_lasso is not None and np.isfinite(z_lasso):
            results["lassosum_z"] = z_lasso
            results["lassosum_pval"] = 2 * norm.sf(abs(z_lasso))
        else:
            results["lassosum_z"] = np.nan
            results["lassosum_pval"] = np.nan
    else:
        results["lassosum_z"] = np.nan
        results["lassosum_pval"] = np.nan

    # ACAT combination
    method_pvals = []
    for key in ["pt_pval", "lassosum_pval"]:
        p = results.get(key)
        if p is not None and np.isfinite(p):
            method_pvals.append(p)

    if method_pvals:
        acat_p = acat_combine(method_pvals)
        if acat_p is not None:
            results["acat_pval"] = acat_p
            # Convert ACAT p-value back to z-score (signed by best method)
            best_z = results.get("pt_z", 0)
            if abs(results.get("lassosum_z", 0) or 0) > abs(best_z or 0):
                best_z = results.get("lassosum_z", 0)
            sign = 1 if (best_z or 0) >= 0 else -1
            results["acat_z"] = sign * norm.isf(acat_p / 2) if acat_p < 1 else 0.0
        else:
            results["acat_pval"] = np.nan
            results["acat_z"] = np.nan
    else:
        results["acat_pval"] = np.nan
        results["acat_z"] = np.nan

    return results


# =============================================================================
# Checkpoint management
# =============================================================================

def save_checkpoint(gwas_name, processed_genes, results_list, out_dir):
    """Save checkpoint with processed genes and partial results."""
    ckpt_path = os.path.join(out_dir, f"checkpoint_{gwas_name}.json")
    ckpt_data = {
        "processed": sorted(processed_genes),
        "n_results": len(results_list),
    }
    with open(ckpt_path, "w") as f:
        json.dump(ckpt_data, f)

    # Also save partial results
    if results_list:
        partial_path = os.path.join(out_dir, f"partial_{gwas_name}.csv")
        pd.DataFrame(results_list).to_csv(partial_path, index=False)


def load_checkpoint(gwas_name, out_dir):
    """Load checkpoint. Returns (set of processed gene IDs, list of partial results)."""
    ckpt_path = os.path.join(out_dir, f"checkpoint_{gwas_name}.json")
    partial_path = os.path.join(out_dir, f"partial_{gwas_name}.csv")

    processed = set()
    results = []

    if os.path.exists(ckpt_path):
        with open(ckpt_path) as f:
            data = json.load(f)
        processed = set(data.get("processed", []))

    if os.path.exists(partial_path):
        results = pd.read_csv(partial_path).to_dict("records")

    return processed, results


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Run TWAS using OTTERS weights against GWAS summary statistics"
    )
    parser.add_argument(
        "--gwas", type=str, default=None,
        help="GWAS name (default: from GWAS_NAME env var). "
             f"Options: {list(GWAS_REGISTRY.keys())}"
    )
    parser.add_argument(
        "--project-root", type=str, default=None,
        help="Project root directory"
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="Resume from checkpoint"
    )
    args = parser.parse_args()

    project_root = args.project_root or get_project_root()
    gwas_name = args.gwas or os.environ.get("GWAS_NAME", "ghodsian")
    gwas_name = gwas_name.lower()

    if gwas_name not in GWAS_REGISTRY:
        print(f"ERROR: Unknown GWAS '{gwas_name}'. Options: {list(GWAS_REGISTRY.keys())}")
        sys.exit(1)

    gwas_info = GWAS_REGISTRY[gwas_name]

    # Paths
    gwas_path = os.path.join(project_root, "GWAS", "MR_Data", "hg19", gwas_info["file"])
    weight_dir = os.path.join(project_root, "data", "broadaway_eqtl", "otters_weights")
    anno_path = os.path.join(project_root, "data", "broadaway_eqtl", "otters_format",
                             "gene_anno.txt")
    ld_dir = os.path.join(project_root, "data", "1kg_eur")
    out_dir = os.path.join(project_root, "RNA-seq", "results", "causal_inference",
                           "otters_broadaway", gwas_name)

    # Validate inputs
    for path, label in [(gwas_path, "GWAS"), (anno_path, "Gene annotation")]:
        if not os.path.exists(path):
            print(f"ERROR: {label} file not found: {path}")
            sys.exit(1)

    os.makedirs(out_dir, exist_ok=True)

    print("=" * 80)
    print(f"OTTERS TWAS — {gwas_name}")
    print("=" * 80)
    print(f"Project root:  {project_root}")
    print(f"GWAS file:     {gwas_path}")
    print(f"GWAS N:        {gwas_info['N']:,}")
    print(f"Weight dir:    {weight_dir}")
    print(f"Gene anno:     {anno_path}")
    print(f"Output dir:    {out_dir}")
    print("=" * 80)

    # Load gene annotation
    print("\nLoading gene annotation...")
    gene_anno = pd.read_csv(anno_path, sep="\t")
    print(f"  {len(gene_anno)} genes across {gene_anno['CHROM'].nunique()} chromosomes")

    # Load GWAS (full — we'll filter by chromosome in the loop)
    print(f"\nLoading GWAS: {gwas_name}...")
    t0 = time.time()
    gwas_full = load_gwas(gwas_path, gwas_info["cols"])
    print(f"  {len(gwas_full):,} GWAS variants loaded ({time.time() - t0:.1f}s)")

    # Resume from checkpoint
    processed_genes = set()
    results_list = []
    if args.resume:
        processed_genes, results_list = load_checkpoint(gwas_name, out_dir)
        if processed_genes:
            print(f"\nResuming: {len(processed_genes)} genes already processed, "
                  f"{len(results_list)} with results")

    # Identify available chromosomes with trained weights
    available_chroms = []
    for c in range(1, 23):
        chr_weight_dir = os.path.join(weight_dir, f"chr{c}")
        if os.path.isdir(chr_weight_dir):
            weight_files = [f for f in os.listdir(chr_weight_dir)
                            if f.endswith(".weights.tsv")]
            if weight_files:
                available_chroms.append(c)

    if not available_chroms:
        print("\nERROR: No trained weights found. Run Script 62 first.")
        sys.exit(1)

    print(f"\nChromosomes with weights: {available_chroms}")

    # Create temp directory
    tmpdir = tempfile.mkdtemp(prefix=f"otters_twas_{gwas_name}_")
    print(f"Temp directory: {tmpdir}")

    t_start = time.time()
    n_success = 0
    n_skipped = 0
    n_no_overlap = 0
    n_error = 0
    gene_counter = 0

    try:
        for chrom in available_chroms:
            chr_weight_dir = os.path.join(weight_dir, f"chr{chrom}")
            pvar_path = os.path.join(ld_dir, f"chr{chrom}_eur.pvar")
            pgen_prefix = os.path.join(ld_dir, f"chr{chrom}_eur")

            if not os.path.exists(pvar_path):
                print(f"\n  WARNING: pvar not found for chr{chrom}, skipping")
                continue

            # Load pvar for this chromosome
            pvar_df = load_pvar_positions(pvar_path)

            # Filter GWAS to this chromosome
            gwas_chr = gwas_full[gwas_full["chr"] == chrom].copy()
            if len(gwas_chr) == 0:
                print(f"\n  WARNING: No GWAS variants on chr{chrom}, skipping")
                continue

            # Get genes on this chromosome
            chr_genes = gene_anno[gene_anno["CHROM"] == chrom]

            # Get weight files for this chromosome
            weight_files = {
                f.replace(".weights.tsv", ""): os.path.join(chr_weight_dir, f)
                for f in os.listdir(chr_weight_dir)
                if f.endswith(".weights.tsv")
            }

            print(f"\n--- Chromosome {chrom}: {len(chr_genes)} genes, "
                  f"{len(weight_files)} weight files, "
                  f"{len(gwas_chr):,} GWAS variants ---")

            for _, gene_row in chr_genes.iterrows():
                gene_id = gene_row["TargetID"]
                gene_name = gene_row["GeneName"]
                gene_counter += 1

                # Skip if already processed
                if gene_id in processed_genes:
                    n_skipped += 1
                    continue

                # Check if weights exist
                if gene_id not in weight_files:
                    processed_genes.add(gene_id)
                    n_no_overlap += 1
                    continue

                # Run TWAS for this gene
                try:
                    result = run_twas_for_gene(
                        weight_files[gene_id], gwas_chr, pvar_df,
                        pgen_prefix, tmpdir
                    )
                except Exception as e:
                    print(f"    ERROR: {gene_id} ({gene_name}): {e}")
                    n_error += 1
                    processed_genes.add(gene_id)
                    continue

                if result is not None:
                    result["gene"] = gene_id
                    result["gene_symbol"] = gene_name
                    result["chr"] = chrom
                    results_list.append(result)
                    n_success += 1
                else:
                    n_no_overlap += 1

                processed_genes.add(gene_id)

                # Progress
                n_proc = len(processed_genes)
                if n_proc % PROGRESS_INTERVAL == 0:
                    elapsed = time.time() - t_start
                    print(f"  [{n_proc} processed] "
                          f"success={n_success}, no_overlap={n_no_overlap}, "
                          f"error={n_error} ({elapsed:.0f}s)")

                # Checkpoint
                if n_proc % CHECKPOINT_INTERVAL == 0:
                    save_checkpoint(gwas_name, processed_genes, results_list, out_dir)

    finally:
        save_checkpoint(gwas_name, processed_genes, results_list, out_dir)
        shutil.rmtree(tmpdir, ignore_errors=True)

    # Build final results
    if not results_list:
        print("\nWARNING: No TWAS results produced. Check weight files and GWAS overlap.")
        sys.exit(0)

    results_df = pd.DataFrame(results_list)

    # Reorder columns
    col_order = [
        "gene", "gene_symbol", "chr", "n_snps",
        "pt_z", "pt_pval",
        "lassosum_z", "lassosum_pval",
        "acat_z", "acat_pval",
    ]
    # Only include columns that exist
    col_order = [c for c in col_order if c in results_df.columns]
    results_df = results_df[col_order]

    # Sort by ACAT p-value
    results_df = results_df.sort_values("acat_pval", ascending=True)

    # FDR correction (Benjamini-Hochberg)
    valid_pvals = results_df["acat_pval"].dropna()
    if len(valid_pvals) > 0:
        from scipy.stats import rankdata
        n_tests = len(valid_pvals)
        ranks = rankdata(valid_pvals)
        fdr = valid_pvals * n_tests / ranks
        # Enforce monotonicity
        fdr_sorted_idx = np.argsort(valid_pvals.values)
        fdr_values = fdr.values.copy()
        fdr_sorted = fdr_values[fdr_sorted_idx]
        for i in range(len(fdr_sorted) - 2, -1, -1):
            fdr_sorted[i] = min(fdr_sorted[i], fdr_sorted[i + 1])
        fdr_values[fdr_sorted_idx] = fdr_sorted
        fdr_values = np.clip(fdr_values, 0, 1)

        results_df["fdr"] = np.nan
        results_df.loc[valid_pvals.index, "fdr"] = fdr_values

    # Save results
    out_path = os.path.join(out_dir, "otters_twas_combined.csv")
    results_df.to_csv(out_path, index=False)

    # Clean up checkpoint files
    for suffix in [f"checkpoint_{gwas_name}.json", f"partial_{gwas_name}.csv"]:
        ckpt_path = os.path.join(out_dir, suffix)
        if os.path.exists(ckpt_path):
            os.remove(ckpt_path)

    # Summary
    t_elapsed = time.time() - t_start
    print("\n" + "=" * 80)
    print(f"SUMMARY — OTTERS TWAS vs {gwas_name}")
    print("=" * 80)
    print(f"  Total genes processed:  {len(processed_genes)}")
    print(f"  TWAS results:           {n_success}")
    print(f"  No GWAS overlap:        {n_no_overlap}")
    print(f"  Errors:                 {n_error}")
    print(f"  Skipped (resumed):      {n_skipped}")

    n_sig_pt = (results_df["pt_pval"] < 0.05).sum() if "pt_pval" in results_df else 0
    n_sig_lasso = (results_df["lassosum_pval"] < 0.05).sum() if "lassosum_pval" in results_df else 0
    n_sig_acat = (results_df["acat_pval"] < 0.05).sum() if "acat_pval" in results_df else 0
    n_fdr = (results_df["fdr"] < 0.05).sum() if "fdr" in results_df else 0

    print(f"\n  Significant (p < 0.05):")
    print(f"    P+T:       {n_sig_pt}")
    print(f"    lassosum:  {n_sig_lasso}")
    print(f"    ACAT:      {n_sig_acat}")
    print(f"    FDR < 0.05:{n_fdr}")

    if n_fdr > 0 and "fdr" in results_df:
        top_hits = results_df[results_df["fdr"] < 0.05].head(20)
        print(f"\n  Top FDR<0.05 hits:")
        for _, row in top_hits.iterrows():
            print(f"    {row['gene_symbol']:<15} z={row.get('acat_z', np.nan):+.3f}  "
                  f"p={row.get('acat_pval', np.nan):.2e}  fdr={row.get('fdr', np.nan):.3e}")

    print(f"\n  Elapsed time:  {t_elapsed:.1f}s ({t_elapsed / 60:.1f} min)")
    print(f"  Output:        {out_path}")

    # Save summary
    summary = {
        "gwas": gwas_name,
        "gwas_n": gwas_info["N"],
        "n_genes_tested": n_success,
        "n_sig_pt_005": int(n_sig_pt),
        "n_sig_lassosum_005": int(n_sig_lasso),
        "n_sig_acat_005": int(n_sig_acat),
        "n_fdr_005": int(n_fdr),
        "elapsed_seconds": round(t_elapsed, 1),
    }
    summary_path = os.path.join(out_dir, f"twas_summary_{gwas_name}.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"  Summary:       {summary_path}")
    print("\nDone.")


if __name__ == "__main__":
    main()
