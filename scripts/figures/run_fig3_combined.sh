#!/bin/bash
#SBATCH --job-name=fig3_combined
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=logs/fig3_combined_%j.out
#SBATCH --error=logs/fig3_combined_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures
mkdir -p logs

echo "=== Fig2 Single-cell & Deconvolution (12 panels) ==="
echo "Start: $(date)"

micromamba run -n rnaseq Rscript fig2_compact.R

echo "=== Done: $(date) ==="
ls -la ../../figures/fig2/fig2_compact.pdf
