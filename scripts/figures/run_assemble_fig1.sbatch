#!/bin/bash
#SBATCH --job-name=assemble_fig1
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/assemble_fig1_%j.out
#SBATCH --error=logs/assemble_fig1_%j.err
set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs figures/main/fig1_atlas_overview

Rscript scripts/figures/assemble_fig1.R
