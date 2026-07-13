#!/bin/bash
#SBATCH --job-name=perm10k
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/perm10k_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/perm10k_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
# Per-gene dual-null convergence permutation FDR (perm10k) on the cleaned-gate
# C2 universe — replaces the stale 2026-05-22 33,943-gene pre-C2 perm10k_fdr.csv.
# Runs 46d with PERM10K_DUAL=TRUE N_PERM=10000; the gated Section 13b block adds
# perm10k_fdr.csv / perm10k_summary.csv / perm10k_top_gene_shuffle.csv. The
# canonical convergence_evidence.csv (L1539, set.seed(42)) is reproduced
# identically since the atlas is unchanged.
set -eo pipefail
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p "$MASLD_PROJECT_ROOT/RNA-seq/logs"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
export PERM10K_DUAL=TRUE
export N_PERM=10000
cd "$MASLD_PROJECT_ROOT/RNA-seq"
echo "[$(date)] START perm10k dual-null (N_PERM=$N_PERM) on C2 universe"
Rscript 46d_convergence_evidence.R
echo "[$(date)] DONE perm10k"
echo "--- perm10k_summary.csv ---"
cat "$MASLD_PROJECT_ROOT/RNA-seq/results/multi_evidence/perm10k_summary.csv"
