#!/bin/bash
#SBATCH --job-name=fig2_regen
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/fig2_regen_%j.out
#SBATCH --error=scripts/figures/logs/fig2_regen_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

source ~/.bashrc
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="${BASE}"

echo "[$(date)] Regenerating Figure 2 (compact) with ashr-shrunk dream results..."
Rscript scripts/figures/fig2_compact.R
echo "[$(date)] Figure 2 DONE"

ls -lh figures/main/fig2_f2_switch/fig2_compact.pdf
