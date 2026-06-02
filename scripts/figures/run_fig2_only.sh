#!/bin/bash
#SBATCH --job-name=liver_fig2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=logs/fig2_regen_%j.log

set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "$MASLD_PROJECT_ROOT"
mkdir -p logs figures

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript scripts/figures/fig2_disease_progression.R
ls -la figures/fig2_disease_progression.pdf
