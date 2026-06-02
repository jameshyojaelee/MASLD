#!/bin/bash
#SBATCH --job-name=mb_bg_cmp
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/56e_cmp_%j.out
#SBATCH --error=logs/56e_cmp_%j.err
set -e

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs results/gwas_atac
Rscript src/56e_compare_backgrounds.R
