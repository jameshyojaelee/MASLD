#!/bin/bash
#SBATCH --job-name=210_progression_coloc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/210_progression_coloc_%j.out
#SBATCH --error=logs/210_progression_coloc_%j.err

# Activate environment BEFORE set -u (activate-binutils references unbound ADDR2LINE)
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -euo pipefail

echo "=== 210_progression_coloc ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"

# Create logs directory if needed
mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

Rscript /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/210_progression_coloc.R

echo "Completed: $(date)"
