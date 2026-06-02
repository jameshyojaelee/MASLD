#!/bin/bash
#SBATCH --partition=bigmem --qos=interactive
#SBATCH --cpus-per-task=16 --mem=200G --time=48:00:00
#SBATCH --job-name=nmf95b_seeds
#SBATCH --output=RNA-seq/logs/nmf95b_seeds_%j.out
#SBATCH --error=RNA-seq/logs/nmf95b_seeds_%j.err
# Bravo Task 8 — 3-seed stability at k=6 AND k=8 (the two rubric-tied k values).
# Seeds {42, 100, 2024} × nrun=50 × {k=6, k=8}  = 6 fresh fits.
# Runs against existing clean cache (mat_nn only) — no matrix rebuild.
# --qos=interactive bypasses the 7TB nslab memory cap.
set -eo pipefail
mkdir -p RNA-seq/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

export NMF_K_LIST="6,8"
export NMF_RUNS=50
# Refit seed 42 at nrun=50 (cache currently holds nrun=20 fit; don't reuse).
export NMF_REUSE_SEED42_CACHE=FALSE
# Canonical clean cache supplies mat_nn (input matrix) only.
export NMF_CACHE_PATH="RNA-seq/results/subtypes/nmf_results_cache_clean.rds"

echo "=== 95b seed stability: k={6,8} × seeds={42,100,2024} × nrun=50 ==="
echo "  Cache in: $NMF_CACHE_PATH"
Rscript RNA-seq/95b_nmf_seed_sensitivity.R
