#!/usr/bin/env bash
#SBATCH --job-name=plot_instaprism
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=cpu
#SBATCH --output=RNA-seq/deconvolution/logs/plot_instaprism_%j.out
#SBATCH --error=RNA-seq/deconvolution/logs/plot_instaprism_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13

# Environment setup
MICROMAMBA=RNA-seq/in-house_MCD_RNAseq/micromamba
ENV_PREFIX=RNA-seq/.mamba/music_deconv
export PYTHONNOUSERSITE=1

echo "=== Generating Creative Plots for InstaPrism (R) ==="
"${MICROMAMBA}" run -p "${ENV_PREFIX}" \
  Rscript RNA-seq/deconvolution/scripts/10_make_creative_plots.R \
  RNA-seq/deconvolution/results \
  RNA-seq/deconvolution/bulk \
  bayesprism

echo ""
echo "=== Generating Extended Analysis Plots for InstaPrism (Python) ==="
"${MICROMAMBA}" run -p "${ENV_PREFIX}" \
  python3 RNA-seq/deconvolution/scripts/extended_deconv_analysis.py \
  --method bayesprism

echo ""
echo "=== Plotting Complete ==="
