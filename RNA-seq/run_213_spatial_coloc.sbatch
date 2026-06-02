#!/bin/bash
#SBATCH --job-name=213_spatial_coloc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/213_spatial_coloc_%j.out
#SBATCH --error=logs/213_spatial_coloc_%j.err

# Activate environment BEFORE set -u (activate-binutils references unbound ADDR2LINE)
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -euo pipefail

echo "=== Script 213: Zonation-stratified COLOC ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo ""

# Create log directory if needed
mkdir -p logs

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

Rscript 213_spatial_coloc.R

echo ""
echo "=== Done: $(date) ==="
