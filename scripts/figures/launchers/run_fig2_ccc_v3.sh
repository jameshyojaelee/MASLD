#!/bin/bash
#SBATCH --job-name=fig2-ccc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=scripts/figures/logs/fig2_ccc_v3_panels_%j.out
#SBATCH --error=scripts/figures/logs/fig2_ccc_v3_panels_%j.err

set -eo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd $MASLD_PROJECT_ROOT
Rscript scripts/figures/fig2_ccc_v3_panels.R
