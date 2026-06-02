#!/bin/bash
#SBATCH --job-name=age_sensitivity
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/age_sensitivity_%j.out
#SBATCH --error=logs/age_sensitivity_%j.err

set -euo pipefail

echo "=== Age Confounding Sensitivity Analysis ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo ""

# Activate environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

# Create log directory if needed
mkdir -p logs

# Run analysis
Rscript 14.3_age_sensitivity.R

echo ""
echo "Completed: $(date)"
