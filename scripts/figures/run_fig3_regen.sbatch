#!/bin/bash
#SBATCH --job-name=fig3_regen
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=scripts/figures/fig3_regen_%j.out
#SBATCH --error=scripts/figures/fig3_regen_%j.err

set -euo pipefail

eval "$(micromamba shell hook --shell=bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "Starting Fig 3 regeneration at $(date)"
Rscript scripts/figures/fig3_progression.R
echo "Completed at $(date)"
ls -la figures/fig3/fig3_progression.pdf
ls -la figures/fig3/panels/panel_*.pdf
