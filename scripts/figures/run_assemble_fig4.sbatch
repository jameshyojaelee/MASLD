#!/bin/bash
#SBATCH --job-name=assemble_fig4
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/assemble_fig4_%j.out
#SBATCH --error=logs/assemble_fig4_%j.err
set -e
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs figures/main/fig4_validation
Rscript scripts/figures/assemble_fig4.R
