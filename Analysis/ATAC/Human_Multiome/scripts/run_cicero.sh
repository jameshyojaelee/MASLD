#!/bin/bash
#SBATCH --job-name=cicero
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=Analysis/ATAC/Human_Multiome/logs/cicero_%j.log
#SBATCH --error=Analysis/ATAC/Human_Multiome/logs/cicero_%j.log

# ==========================================================================
# Cicero co-accessibility analysis pipeline
#
# Step 1 (Python): Extract peak matrix from SnapATAC2 h5ad, export for R
# Step 2 (R):      Run Cicero co-accessibility + differential + plots
# ==========================================================================

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/ATAC/Human_Multiome/logs
mkdir -p Analysis/ATAC/Human_Multiome/results/cicero

echo "============================================================"
echo "Cicero Co-accessibility Pipeline"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Date:   $(date)"
echo "============================================================"

# --- Step 1: Prepare data (Python) ----------------------------------------
echo ""
echo ">>> Step 1: Preparing peak matrix (Python)..."
echo ""

eval "$(micromamba shell hook -s bash)"
micromamba activate snapatac2

python Analysis/ATAC/Human_Multiome/scripts/10_cicero_prepare.py

echo ""
echo ">>> Step 1 complete."
echo ""

# --- Step 2: Run Cicero (R) -----------------------------------------------
echo ">>> Step 2: Running Cicero analysis (R)..."
echo ""

micromamba activate rnaseq

Rscript Analysis/ATAC/Human_Multiome/scripts/11_cicero_analysis.R

# --- Step 3: Differential co-accessibility (R) ----------------------------
echo ">>> Step 3: Differential co-accessibility..."
echo ""

Rscript Analysis/ATAC/Human_Multiome/scripts/12_cicero_differential.R

echo ""
echo "============================================================"
echo "Pipeline complete."
echo "Results: Analysis/ATAC/Human_Multiome/results/cicero/"
echo "Figures: figures/supplementary/epigenomic/"
echo "============================================================"
