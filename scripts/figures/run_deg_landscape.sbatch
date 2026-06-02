#!/bin/bash
#SBATCH --job-name=deg_landscape
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/deg_landscape_%j.log

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
source ~/.bashrc
micromamba activate rnaseq

Rscript scripts/figures/figS_deg_threshold_landscape.R
