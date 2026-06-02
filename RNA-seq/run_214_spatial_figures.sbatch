#!/bin/bash
#SBATCH --job-name=214_spatial_figures
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/214_spatial_figures_%j.out
#SBATCH --error=logs/214_spatial_figures_%j.err

# Activate environment BEFORE set -u (activate-binutils references unbound ADDR2LINE)
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -euo pipefail

echo "=== Script 214: Spatial-stratified COLOC Figures ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo ""

# Create log directory if needed
mkdir -p logs

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

Rscript 214_spatial_figures.R

echo ""
echo "=== Done: $(date) ==="
