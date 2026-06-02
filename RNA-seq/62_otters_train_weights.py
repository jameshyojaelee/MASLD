#!/usr/bin/env python3
"""
62_otters_train_weights.py
Phase 2A-Step2: Train TWAS prediction weights from Broadaway liver eQTL summary
statistics + 1KG EUR LD reference.

Since OTTERS is not available as a pip package, this script implements the core
weight-training methods directly:
  - P+T (Pruning + Thresholding): Select SNPs by eQTL p-value, LD-prune (r² < 0.1),
    use marginal effect sizes as weights.
  - lassosum: Coordinate descent on the penalised regression objective using summary
    statistics and an external LD matrix (Mak et al., 2017).

Input (from Script 61):
  - data/broadaway_eqtl/otters_format/chr{1-22}_broadaway.txt  (OTTERS sseQTL)
  - data/broadaway_eqtl/otters_format/gene_anno.txt            (gene annotation)
  - data/1kg_eur/chr{1-22}_eur.{pgen,psam,pvar}                (1KG EUR Phase 1, hg19)

Output:
  - data/broadaway_eqtl/otters_weights/chr{CHR}/  — per-gene weight files
    Format: SNP_ID, POS, A1, A2, Weight_PT, Weight_lassosum

Usage:
  python 62_otters_train_weights.py --chr 21
  # Or via SLURM array:
  sbatch run_otters_train.sbatch   # --array=1-22

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
from scipy.stats import norm
from scipy.linalg import cho_factor, cho_solve


# =============================================================================
# Configuration
# =============================================================================

CIS_WINDOW = 1_000_000  # ±1 Mb from gene body
PT_PVAL_THRESHOLD = 0.01  # p-value threshold for P+T
PT_R2_THRESHOLD = 0.1  # LD r² pruning threshold for P+T
LASSOSUM_S_VALUES = [0.2, 0.5, 0.9]  # lassosum shrinkage grid
LASSOSUM_LAMBDA_VALUES = [0.001, 0.005, 0.01, 0.05]  # lassosum penalty grid
LASSOSUM_MAX_ITER = 500  # max iterations for coordinate descent
LASSOSUM_TOL = 1e-6  # convergence tolerance
MIN_SNPS_PER_GENE = 5  # minimum SNPs to attempt training
CHECKPOINT_INTERVAL = 500  # save checkpoint every N genes
PROGRESS_INTERVAL = 100  # print progress every N genes
N_EQTL = 1183  # Broadaway sample size


def get_project_root():
    """Get project root from env or default."""
    return os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
    )


# =============================================================================
# LD computation via plink2
# =============================================================================

def find_plink2():
    """Find plink2 binary — try PATH first, then module load."""
    # Check if already in PATH
    result = shutil.which("plink2")
    if result:
        return result
    # Try module load
    try:
        subprocess.run(
            "module load PLINK/2.00a3.7-gfbf-2023a && which plink2",
            shell=True, capture_output=True, text=True, check=True,
            executable="/bin/bash"
        )
        # If module load works, we'll use it in our commands
        return "plink2"
    except subprocess.CalledProcessError:
        pass
    raise RuntimeError(
        "plink2 not found. Ensure 'module load PLINK/2.00a3.7-gfbf-2023a' works."
    )


def compute_ld_matrix(pgen_prefix, snp_ids, tmpdir):
    """Compute LD r matrix for a set of SNPs using plink2.

    Args:
        pgen_prefix: Path prefix for pgen/psam/pvar files (without extension)
        snp_ids: List of SNP IDs (matching pvar ID column, e.g., '21:9411618:C:T')
        tmpdir: Temporary directory for intermediate files

    Returns:
        ld_matrix: np.ndarray of shape (n_snps, n_snps) — Pearson r (not r²)
        kept_ids: List of SNP IDs in the matrix (some may be dropped by plink2)
    """
    if len(snp_ids) < 2:
        return np.array([[1.0]]), snp_ids

    # Write SNP list
    snp_file = os.path.join(tmpdir, "snps.txt")
    with open(snp_file, "w") as f:
        for sid in snp_ids:
            f.write(sid + "\n")

    out_prefix = os.path.join(tmpdir, "ld_out")

    # Run plink2 --r-unphased square (produces .unphased.vcor matrix)
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

    # Check for output file — plink2 may produce .unphased.vcor or .unphased.vcor1
    vcor_file = None
    for suffix in [".unphased.vcor", ".unphased.vcor1"]:
        candidate = out_prefix + suffix
        if os.path.exists(candidate):
            vcor_file = candidate
            break

    if vcor_file is None:
        # Fallback: try --r square (older plink2 versions)
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
            f"plink2 LD computation failed. stderr:\n{result.stderr}\n"
            f"stdout:\n{result.stdout}"
        )

    # Read LD matrix (tab-separated square matrix, no header)
    ld_matrix = np.loadtxt(vcor_file, delimiter="\t")

    # Read which SNPs were kept (from .vars file if present, else from the ID file)
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
        # Assume all SNPs were kept if no .vars file
        kept_ids = list(snp_ids)
        # Sanity check dimension
        if ld_matrix.shape[0] != len(kept_ids):
            raise RuntimeError(
                f"LD matrix dim {ld_matrix.shape} != SNP count {len(kept_ids)}"
            )

    return ld_matrix, kept_ids


# =============================================================================
# P+T (Pruning + Thresholding) weight training
# =============================================================================

def train_pt_weights(zscore, ld_r, pval_threshold=PT_PVAL_THRESHOLD,
                     r2_threshold=PT_R2_THRESHOLD):
    """Train P+T weights from summary statistics.

    Steps:
      1. Convert Z-scores to two-sided p-values.
      2. Filter SNPs with p < threshold.
      3. LD-clump: greedily select SNPs by ascending p-value, removing any
         remaining SNP with r² > threshold relative to the selected SNP.
      4. Return marginal effect sizes (Z / sqrt(N)) as weights for selected SNPs.

    Args:
        zscore: np.ndarray of Z-scores (n_snps,)
        ld_r: np.ndarray LD correlation matrix (n_snps, n_snps)
        pval_threshold: p-value threshold for inclusion
        r2_threshold: LD r² threshold for pruning

    Returns:
        weights: np.ndarray of weights (n_snps,), zero for pruned SNPs
    """
    n = len(zscore)
    pvals = 2 * norm.sf(np.abs(zscore))
    weights = np.zeros(n)

    # Indices passing p-value threshold, sorted by ascending p-value
    pass_mask = pvals < pval_threshold
    if not np.any(pass_mask):
        # Relax threshold: take top 10 SNPs by p-value
        n_take = min(10, n)
        top_idx = np.argsort(pvals)[:n_take]
        pass_mask = np.zeros(n, dtype=bool)
        pass_mask[top_idx] = True

    candidates = np.where(pass_mask)[0]
    candidates = candidates[np.argsort(pvals[candidates])]

    # Greedy LD clumping
    selected = []
    removed = set()
    ld_r2 = ld_r ** 2

    for idx in candidates:
        if idx in removed:
            continue
        selected.append(idx)
        # Remove all candidates in LD with this SNP
        for other in candidates:
            if other != idx and other not in removed:
                if ld_r2[idx, other] > r2_threshold:
                    removed.add(other)

    # Assign weights: marginal beta ~ Z / sqrt(N)
    for idx in selected:
        weights[idx] = zscore[idx] / np.sqrt(N_EQTL)

    return weights


# =============================================================================
# lassosum weight training
# =============================================================================

def train_lassosum_weights(zscore, ld_r, s_values=LASSOSUM_S_VALUES,
                           lambda_values=LASSOSUM_LAMBDA_VALUES,
                           max_iter=LASSOSUM_MAX_ITER, tol=LASSOSUM_TOL):
    """Train lassosum weights from summary statistics + LD.

    Implements the lassosum coordinate descent (Mak et al., 2017):
      Minimise: (1-s)*beta'*R*beta - 2*beta'*r_zy + lambda*||beta||_1
    where r_zy = Z / sqrt(N) are the marginal correlations.

    Uses a grid of (s, lambda) and selects the model with best pseudo-R²
    (approximate in-sample prediction).

    Args:
        zscore: np.ndarray of Z-scores (n_snps,)
        ld_r: np.ndarray LD correlation matrix (n_snps, n_snps)
        s_values: shrinkage parameter grid (0 < s < 1)
        lambda_values: L1 penalty grid
        max_iter: maximum iterations
        tol: convergence tolerance

    Returns:
        weights: np.ndarray of weights (n_snps,)
    """
    n_snps = len(zscore)
    r_zy = zscore / np.sqrt(N_EQTL)  # marginal correlations

    best_weights = np.zeros(n_snps)
    best_pseudo_r2 = -np.inf

    for s in s_values:
        # Regularised LD: (1 - s) * R + s * I
        R_reg = (1 - s) * ld_r + s * np.eye(n_snps)

        for lam in lambda_values:
            beta = _coordinate_descent_lasso(R_reg, r_zy, lam, max_iter, tol)

            # Pseudo R² = 2 * beta' * r_zy - beta' * R_reg * beta
            pseudo_r2 = 2 * np.dot(beta, r_zy) - np.dot(beta, R_reg @ beta)

            if pseudo_r2 > best_pseudo_r2:
                best_pseudo_r2 = pseudo_r2
                best_weights = beta.copy()

    return best_weights


def _coordinate_descent_lasso(R, r_zy, lam, max_iter, tol):
    """Coordinate descent for L1-penalised regression with summary stats.

    Solves: min_beta  beta' R beta - 2 beta' r_zy + lam * ||beta||_1
    """
    n = len(r_zy)
    beta = np.zeros(n)
    diag_R = np.diag(R).copy()
    # Avoid division by zero
    diag_R[diag_R < 1e-10] = 1e-10

    for iteration in range(max_iter):
        beta_old = beta.copy()
        for j in range(n):
            # Partial residual
            r_j = r_zy[j] - np.dot(R[j, :], beta) + R[j, j] * beta[j]
            # Soft-thresholding
            beta[j] = _soft_threshold(r_j, lam) / diag_R[j]

        # Check convergence
        if np.max(np.abs(beta - beta_old)) < tol:
            break

    return beta


def _soft_threshold(x, lam):
    """Soft-thresholding operator for L1 penalty."""
    if x > lam:
        return x - lam
    elif x < -lam:
        return x + lam
    else:
        return 0.0


# =============================================================================
# Gene processing
# =============================================================================

def load_eqtl_chr(eqtl_path):
    """Load OTTERS-formatted eQTL summary stats for one chromosome.

    Returns DataFrame with columns: CHROM, POS, A1, A2, Zscore, TargetID, N
    """
    df = pd.read_csv(eqtl_path, sep="\t")
    df["POS"] = df["POS"].astype(np.int64)
    df["Zscore"] = df["Zscore"].astype(np.float64)
    return df


def load_gene_anno(anno_path, chrom):
    """Load gene annotation and filter to the specified chromosome."""
    df = pd.read_csv(anno_path, sep="\t")
    df = df[df["CHROM"] == chrom].copy()
    df["GeneStart"] = df["GeneStart"].astype(np.int64)
    df["GeneEnd"] = df["GeneEnd"].astype(np.int64)
    return df


def load_pvar_positions(pvar_path):
    """Load pvar file to get variant ID -> position mapping.

    Returns DataFrame with columns: POS, ID, REF, ALT
    """
    df = pd.read_csv(
        pvar_path, sep="\t", comment="#", header=None,
        names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"],
        usecols=[1, 2, 3, 4],
    )
    df.columns = ["POS", "ID", "REF", "ALT"]
    df["POS"] = df["POS"].astype(np.int64)
    # Biallelic SNPs only
    snp_mask = (df["REF"].str.len() == 1) & (df["ALT"].str.len() == 1)
    df = df[snp_mask].copy()
    return df


def process_gene(gene_id, gene_name, gene_start, gene_end, eqtl_df,
                 pvar_df, pgen_prefix, chrom, tmpdir):
    """Process a single gene: extract cis-window SNPs, compute LD, train weights.

    Returns:
        result_df: DataFrame with columns [SNP_ID, POS, A1, A2, Weight_PT, Weight_lassosum]
                   or None if the gene is skipped.
    """
    # Define cis-window
    window_start = max(0, gene_start - CIS_WINDOW)
    window_end = gene_end + CIS_WINDOW

    # Extract eQTL summary stats in the cis-window for this gene
    mask = (
        (eqtl_df["TargetID"] == gene_id) &
        (eqtl_df["POS"] >= window_start) &
        (eqtl_df["POS"] <= window_end)
    )
    gene_eqtl = eqtl_df[mask].copy()

    if len(gene_eqtl) < MIN_SNPS_PER_GENE:
        return None

    # Match eQTL SNPs to pvar IDs (by position + alleles)
    # eQTL has A1 (effect/ALT in 1KG frame), A2 (other/REF in 1KG frame)
    # pvar has REF, ALT
    merged = gene_eqtl.merge(pvar_df, on="POS", how="inner", suffixes=("", "_pvar"))

    # Direct match: eQTL A1 == pvar ALT and eQTL A2 == pvar REF (no flip)
    direct = (merged["A1"] == merged["ALT"]) & (merged["A2"] == merged["REF"])
    # Flipped: eQTL A1 == pvar REF and eQTL A2 == pvar ALT (need Z-score flip)
    flipped = (merged["A1"] == merged["REF"]) & (merged["A2"] == merged["ALT"])

    matched = merged[direct | flipped].copy()
    if len(matched) < MIN_SNPS_PER_GENE:
        return None

    # Flip Z-scores where alleles are reversed
    flip_mask = flipped[direct | flipped].values
    matched.loc[flip_mask, "Zscore"] = -matched.loc[flip_mask, "Zscore"]

    # Deduplicate by position
    matched = matched.drop_duplicates(subset=["POS"], keep="first")

    if len(matched) < MIN_SNPS_PER_GENE:
        return None

    # Get SNP IDs for plink2
    snp_ids = matched["ID"].tolist()
    zscores = matched["Zscore"].values.astype(np.float64)

    # Compute LD matrix via plink2
    try:
        ld_r, kept_ids = compute_ld_matrix(pgen_prefix, snp_ids, tmpdir)
    except RuntimeError as e:
        print(f"    WARNING: LD computation failed for {gene_id}: {e}")
        return None

    # Subset to kept SNPs (plink2 may drop some)
    if len(kept_ids) != len(snp_ids):
        kept_set = set(kept_ids)
        keep_mask = matched["ID"].isin(kept_set)
        matched = matched[keep_mask].copy()
        zscores = matched["Zscore"].values.astype(np.float64)

    if len(matched) < MIN_SNPS_PER_GENE:
        return None

    n_snps = len(matched)

    # Sanity check LD matrix dimensions
    if ld_r.shape[0] != n_snps or ld_r.shape[1] != n_snps:
        print(f"    WARNING: LD dim mismatch for {gene_id}: "
              f"matrix {ld_r.shape} vs {n_snps} SNPs. Skipping.")
        return None

    # Fix any NaN/inf in LD matrix
    if np.any(~np.isfinite(ld_r)):
        ld_r = np.nan_to_num(ld_r, nan=0.0, posinf=1.0, neginf=-1.0)
    # Ensure diagonal is 1.0
    np.fill_diagonal(ld_r, 1.0)

    # Train P+T weights
    try:
        w_pt = train_pt_weights(zscores, ld_r)
    except Exception as e:
        print(f"    WARNING: P+T failed for {gene_id}: {e}")
        w_pt = np.zeros(n_snps)

    # Train lassosum weights
    try:
        w_lasso = train_lassosum_weights(zscores, ld_r)
    except Exception as e:
        print(f"    WARNING: lassosum failed for {gene_id}: {e}")
        w_lasso = np.zeros(n_snps)

    # Build result DataFrame
    result_df = pd.DataFrame({
        "SNP_ID": matched["ID"].values,
        "POS": matched["POS"].values,
        "A1": matched["A1"].values,
        "A2": matched["A2"].values,
        "Weight_PT": w_pt,
        "Weight_lassosum": w_lasso,
    })

    return result_df


# =============================================================================
# Checkpoint management
# =============================================================================

def save_checkpoint(chrom, processed_genes, out_dir):
    """Save checkpoint: list of processed gene IDs."""
    ckpt_path = os.path.join(out_dir, f"checkpoint_chr{chrom}.json")
    with open(ckpt_path, "w") as f:
        json.dump({"processed": sorted(processed_genes)}, f)


def load_checkpoint(chrom, out_dir):
    """Load checkpoint if it exists. Returns set of processed gene IDs."""
    ckpt_path = os.path.join(out_dir, f"checkpoint_chr{chrom}.json")
    if os.path.exists(ckpt_path):
        with open(ckpt_path) as f:
            data = json.load(f)
        return set(data.get("processed", []))
    return set()


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Train TWAS prediction weights from Broadaway eQTL + 1KG EUR LD"
    )
    parser.add_argument(
        "--chr", type=int, required=True,
        help="Chromosome to process (1-22)"
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

    chrom = args.chr
    project_root = args.project_root or get_project_root()

    # Paths
    eqtl_dir = os.path.join(project_root, "data", "broadaway_eqtl", "otters_format")
    ld_dir = os.path.join(project_root, "data", "1kg_eur")
    out_dir = os.path.join(project_root, "data", "broadaway_eqtl", "otters_weights",
                           f"chr{chrom}")

    eqtl_path = os.path.join(eqtl_dir, f"chr{chrom}_broadaway.txt")
    anno_path = os.path.join(eqtl_dir, "gene_anno.txt")
    pvar_path = os.path.join(ld_dir, f"chr{chrom}_eur.pvar")
    pgen_prefix = os.path.join(ld_dir, f"chr{chrom}_eur")

    # Validate inputs
    for path, label in [(eqtl_path, "eQTL"), (anno_path, "annotation"),
                        (pvar_path, "pvar")]:
        if not os.path.exists(path):
            print(f"ERROR: {label} file not found: {path}")
            sys.exit(1)

    # Create output directory
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 80)
    print(f"OTTERS WEIGHT TRAINING — Chromosome {chrom}")
    print("=" * 80)
    print(f"Project root:  {project_root}")
    print(f"eQTL file:     {eqtl_path}")
    print(f"Gene anno:     {anno_path}")
    print(f"LD prefix:     {pgen_prefix}")
    print(f"Output dir:    {out_dir}")
    print(f"P+T threshold: p < {PT_PVAL_THRESHOLD}, r² < {PT_R2_THRESHOLD}")
    print(f"lassosum grid: s={LASSOSUM_S_VALUES}, lambda={LASSOSUM_LAMBDA_VALUES}")
    print(f"Cis-window:    ±{CIS_WINDOW / 1e6:.1f} Mb")
    print(f"Min SNPs:      {MIN_SNPS_PER_GENE}")
    print("=" * 80)

    # Verify plink2 is accessible
    try:
        find_plink2()
        print("plink2: OK")
    except RuntimeError as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    # Load data
    print(f"\nLoading eQTL data for chr{chrom}...")
    t0 = time.time()
    eqtl_df = load_eqtl_chr(eqtl_path)
    print(f"  {len(eqtl_df):,} eQTL records, "
          f"{eqtl_df['TargetID'].nunique()} genes ({time.time() - t0:.1f}s)")

    print(f"Loading gene annotation for chr{chrom}...")
    gene_anno = load_gene_anno(anno_path, chrom)
    print(f"  {len(gene_anno)} genes on chr{chrom}")

    print(f"Loading pvar for chr{chrom}...")
    pvar_df = load_pvar_positions(pvar_path)
    print(f"  {len(pvar_df):,} SNPs in 1KG")

    # Resume from checkpoint if requested
    processed_genes = set()
    if args.resume:
        processed_genes = load_checkpoint(chrom, out_dir)
        if processed_genes:
            print(f"\nResuming: {len(processed_genes)} genes already processed")

    # Process genes
    total_genes = len(gene_anno)
    n_success = 0
    n_skipped = 0
    n_too_few = 0
    n_error = 0
    t_start = time.time()

    # Create a persistent temp directory for plink2 intermediate files
    tmpdir = tempfile.mkdtemp(prefix=f"otters_chr{chrom}_")
    print(f"\nTemp directory: {tmpdir}")
    print(f"\nProcessing {total_genes} genes on chr{chrom}...\n")

    try:
        for i, (_, gene_row) in enumerate(gene_anno.iterrows()):
            gene_id = gene_row["TargetID"]
            gene_name = gene_row["GeneName"]
            gene_start = gene_row["GeneStart"]
            gene_end = gene_row["GeneEnd"]

            # Skip if already processed (checkpoint)
            if gene_id in processed_genes:
                n_skipped += 1
                continue

            # Clean tmpdir for each gene to avoid file accumulation
            for f in os.listdir(tmpdir):
                fpath = os.path.join(tmpdir, f)
                try:
                    os.remove(fpath)
                except OSError:
                    pass

            # Process gene
            try:
                result = process_gene(
                    gene_id, gene_name, gene_start, gene_end,
                    eqtl_df, pvar_df, pgen_prefix, chrom, tmpdir
                )
            except Exception as e:
                print(f"  ERROR processing {gene_id} ({gene_name}): {e}")
                n_error += 1
                processed_genes.add(gene_id)
                continue

            if result is None:
                n_too_few += 1
            else:
                # Save gene weights
                gene_out_path = os.path.join(out_dir, f"{gene_id}.weights.tsv")
                result.to_csv(gene_out_path, sep="\t", index=False)
                n_success += 1

            processed_genes.add(gene_id)

            # Progress report
            n_done = i + 1
            if n_done % PROGRESS_INTERVAL == 0 or n_done == total_genes:
                elapsed = time.time() - t_start
                rate = n_done / elapsed if elapsed > 0 else 0
                eta = (total_genes - n_done) / rate if rate > 0 else 0
                print(f"  [{n_done}/{total_genes}] "
                      f"success={n_success}, too_few={n_too_few}, error={n_error} "
                      f"({elapsed:.0f}s elapsed, ~{eta:.0f}s remaining)")

            # Checkpoint
            if n_done % CHECKPOINT_INTERVAL == 0:
                save_checkpoint(chrom, processed_genes, out_dir)

    finally:
        # Always save final checkpoint and clean up
        save_checkpoint(chrom, processed_genes, out_dir)
        shutil.rmtree(tmpdir, ignore_errors=True)

    # Summary
    t_elapsed = time.time() - t_start
    print("\n" + "=" * 80)
    print(f"SUMMARY — Chromosome {chrom}")
    print("=" * 80)
    print(f"  Total genes:        {total_genes}")
    print(f"  Successful:         {n_success}")
    print(f"  Too few SNPs:       {n_too_few}")
    print(f"  Errors:             {n_error}")
    print(f"  Skipped (resumed):  {n_skipped}")
    print(f"  Elapsed time:       {t_elapsed:.1f}s ({t_elapsed / 60:.1f} min)")

    # Write summary file
    summary = {
        "chrom": chrom,
        "total_genes": total_genes,
        "n_success": n_success,
        "n_too_few_snps": n_too_few,
        "n_error": n_error,
        "n_skipped_resume": n_skipped,
        "elapsed_seconds": round(t_elapsed, 1),
        "parameters": {
            "cis_window": CIS_WINDOW,
            "pt_pval_threshold": PT_PVAL_THRESHOLD,
            "pt_r2_threshold": PT_R2_THRESHOLD,
            "lassosum_s_values": LASSOSUM_S_VALUES,
            "lassosum_lambda_values": LASSOSUM_LAMBDA_VALUES,
            "min_snps_per_gene": MIN_SNPS_PER_GENE,
        }
    }
    summary_path = os.path.join(out_dir, f"training_summary_chr{chrom}.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n  Summary written to: {summary_path}")

    # Count weight files
    weight_files = [f for f in os.listdir(out_dir) if f.endswith(".weights.tsv")]
    print(f"  Weight files: {len(weight_files)}")
    print("\nDone.")


if __name__ == "__main__":
    main()
