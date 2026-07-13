#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --job-name=loo_post
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/loo_cv_post_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/loo_cv_post_%j.err

# Run AFTER all 8 dream LOO-CV jobs complete.
# Aggregates LOO results + computes cross-study AUROC + generates figures.

set -euo pipefail

export PATH="/gpfs/commons/home/jameslee/miniforge3/condabin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "=== Step 1: Aggregating LOO-CV results ==="
Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/aggregate_loo_cv.R

echo ""
echo "=== Step 2: Cross-study AUROC prediction ==="
Rscript scripts/cross_study_auroc.R

echo ""
echo "=== All post-LOO-CV processing complete ==="
