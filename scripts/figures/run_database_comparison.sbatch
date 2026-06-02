#!/bin/bash
#SBATCH --job-name=figS_db_comparison
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/figS_db_comparison_%j.out
#SBATCH --error=logs/figS_db_comparison_%j.err

set -euo pipefail

echo "=== figS_database_comparison ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Date: $(date)"

# Activate environment
eval "$(micromamba shell hook --shell=bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

mkdir -p "${MASLD_PROJECT_ROOT}/scripts/figures/logs"

Rscript "${MASLD_PROJECT_ROOT}/scripts/figures/figS_database_comparison.R"

echo "=== DONE: $(date) ==="
