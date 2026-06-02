#!/bin/bash
#SBATCH --job-name=variance_partition
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=logs/variance_partition_%j.out
#SBATCH --error=logs/variance_partition_%j.err

set -euo pipefail

echo "=== Variance Partition Analysis ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Start: $(date)"

# Create log directory if needed
mkdir -p logs

# Activate environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

Rscript "${MASLD_PROJECT_ROOT}/scripts/figures/figS_variance_partition.R"

echo "End: $(date)"
