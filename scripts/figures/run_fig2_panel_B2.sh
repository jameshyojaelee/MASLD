#!/bin/bash
#SBATCH --job-name=fig2B2
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=logs/fig2B2_%j.log
#SBATCH --error=logs/fig2B2_%j.err

set -eo pipefail
mkdir -p logs

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript scripts/figures/fig2_panel_B2.R
