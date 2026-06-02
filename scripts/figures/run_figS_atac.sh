#!/bin/bash
#SBATCH --job-name=figS_atac
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=scripts/figures/logs/figS_atac_%j.log
#SBATCH --error=scripts/figures/logs/figS_atac_%j.log

# ==========================================================================
# Supplementary ATAC-seq figure panels
# Generates 16 individual panels in figures/supplementary/figS05_epigenomic_spatial/
#
# Panels 01-08: Python (h5ad-dependent: QC, UMAPs, dot plot, label transfer)
# Panels 09-16: R (CSV-dependent: volcanos, heatmaps, correlations)
# ==========================================================================

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p scripts/figures/logs
mkdir -p figures/supplementary/figS05_epigenomic_spatial

echo "============================================================"
echo "Supplementary ATAC Figure Panels"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Date:   $(date)"
echo "============================================================"

# --- Python panels (01-08) ------------------------------------------------
echo ""
echo ">>> Running Python panels (01-08)..."
echo ""

# Use snapatac2 env for h5ad compatibility
eval "$(micromamba shell hook -s bash)"
micromamba activate snapatac2

python scripts/figures/figS_atac_python_panels.py

echo ""
echo ">>> Python panels done."
echo ""

# --- R panels (09-16) ----------------------------------------------------
echo ">>> Running R panels (09-16)..."
echo ""

micromamba activate rnaseq

Rscript scripts/figures/figS_atac.R

echo ""
echo "============================================================"
echo "All 16 panels complete."
echo "Output: figures/supplementary/figS05_epigenomic_spatial/"
echo "============================================================"
echo ""
ls -la figures/supplementary/figS05_epigenomic_spatial/
