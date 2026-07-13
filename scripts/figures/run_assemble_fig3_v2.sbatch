#!/bin/bash
#SBATCH --job-name=assemble_fig3
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/assemble_fig3_%j.out
#SBATCH --error=logs/assemble_fig3_%j.err
set -e
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs figures/main/fig2_genetics
Rscript scripts/figures/assemble_fig3_v2.R
