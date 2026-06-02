#!/bin/bash
#SBATCH --job-name=hardening_e4_pca
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --qos=interactive
#SBATCH --output=scripts/logs/hardening_e4_pca_%j.out
#SBATCH --error=scripts/logs/hardening_e4_pca_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "=== Task 1: E4 Drug Tier Verification ==="
echo "Start: $(date)"
micromamba run -n rnaseq Rscript RNA-seq/results/drug_repurposing/e4_tier_verification.R
echo "Task 1 done: $(date)"

echo ""
echo "=== Task 2: Batch-effect PCA Visualization ==="
echo "Start: $(date)"
micromamba run -n rnaseq Rscript scripts/figures/figS_batch_pca.R
echo "Task 2 done: $(date)"

echo ""
echo "=== All tasks complete ==="
echo "End: $(date)"
