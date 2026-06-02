#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --job-name=figS_coloc_abf
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/figS_coloc_abf_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/figS_coloc_abf_%j.err
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

mkdir -p scripts/figures/logs

eval "$(micromamba shell hook --shell=bash)"
micromamba activate finemapping

Rscript scripts/figures/figS04_coloc.R
