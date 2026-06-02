#!/bin/bash
#SBATCH --job-name=fig1_regen
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/fig1_regen_%j.out
#SBATCH --error=scripts/figures/logs/fig1_regen_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

source ~/.bashrc
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="${BASE}"

echo "[$(date)] Regenerating Figure 1 (compact) with ashr-shrunk dream results..."
Rscript scripts/figures/fig1_compact.R
echo "[$(date)] Figure 1 DONE"

ls -lh figures/main/fig1_atlas_overview/fig1_compact.pdf
