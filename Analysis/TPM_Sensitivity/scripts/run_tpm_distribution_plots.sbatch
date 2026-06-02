#!/usr/bin/env bash
#SBATCH --job-name=tpm_dist_plots
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --partition=cpu
#SBATCH --output=RNA-seq/TPM_analysis/logs/tpm_dist_plots_%j.out
#SBATCH --error=RNA-seq/TPM_analysis/logs/tpm_dist_plots_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# Create log directory
mkdir -p RNA-seq/TPM_analysis/logs

# Environment setup
MICROMAMBA=RNA-seq/in-house_MCD_RNAseq/micromamba
ENV_PREFIX=RNA-seq/.mamba/music_deconv
export PYTHONNOUSERSITE=1
export MPLBACKEND=Agg

echo "=== Generating Publication-Quality TPM Distribution Plots ==="
echo "Start time: $(date)"
echo ""

"${MICROMAMBA}" run -p "${ENV_PREFIX}" \
  Rscript RNA-seq/scripts/plot_tpm_distributions_publication.R

echo ""
echo "=== Plotting Complete ==="
echo "End time: $(date)"
echo ""
echo "Output files:"
ls -la RNA-seq/plots/tpm_distribution*.pdf 2>/dev/null || echo "No PDF files found"
