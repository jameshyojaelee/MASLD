#!/bin/bash
#SBATCH --job-name=fm_figures
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/08_finemapping_figures_%j.out
#SBATCH --error=logs/08_finemapping_figures_%j.err

# =============================================================================
# 08_finemapping_figures.sh
# Generate publication-quality fine-mapping summary figures
# =============================================================================

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs

echo "=== Fine-mapping figures ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID:-local}"

# Activate finemapping conda env
eval "$(micromamba shell hook -s bash)"
micromamba activate finemapping

echo "R version: $(R --version | head -1)"

Rscript src/08_finemapping_figures.R

echo "=== Done: $(date) ==="
