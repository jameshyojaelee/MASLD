#!/bin/bash
#SBATCH --job-name=fisher
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/46e_fisher_%j.out
#SBATCH --error=logs/46e_fisher_%j.err

set -euo pipefail
mkdir -p logs

echo "=== 46e Fisher Combined Test + COLOC Enzyme Sensitivity ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"

micromamba run -n rnaseq Rscript 46e_fisher_combined_test.R

echo "End: $(date)"
