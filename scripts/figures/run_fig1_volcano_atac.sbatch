#!/bin/bash
#SBATCH --job-name=fig1_volcano_atac
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/fig1_volcano_atac_%j.out
#SBATCH --error=logs/fig1_volcano_atac_%j.err
set -e
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs figures/main/fig1_atlas_overview/panels
Rscript scripts/figures/fig1_volcano_atac_overlay.R
