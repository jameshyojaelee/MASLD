#!/usr/bin/env bash
#SBATCH --job-name=sc_figures
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=scripts/figures/logs/sc_figures_%j.out
#SBATCH --error=scripts/figures/logs/sc_figures_%j.err

set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "=== Generating Figure 3 (deconvolution + sc panels e-h) ==="
Rscript scripts/figures/fig3_deconvolution_hero.R

echo ""
echo "=== Generating Figure 5 (single-cell atlas) ==="
Rscript scripts/figures/fig5_singlecell.R

echo ""
echo "=== Done ==="
ls -la scripts/figures/output/fig3_deconvolution_hero.pdf
ls -la scripts/figures/output/fig5_singlecell.pdf
echo "End: $(date)"
