#!/bin/bash
#SBATCH --job-name=fig3_panelJ
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/fig3_panelJ_%j.out
#SBATCH --error=logs/fig3_panelJ_%j.err
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs figures/main/fig3_regulatory_architecture/panels
Rscript scripts/figures/fig3_panel_J_disease_regulons_celltype.R
