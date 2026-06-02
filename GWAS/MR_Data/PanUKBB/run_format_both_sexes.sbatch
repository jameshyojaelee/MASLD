#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --job-name=B1_BS_fmt
#SBATCH --output=GWAS/MR_Data/PanUKBB/sex_stratified/logs/format_BS_%j.out
#SBATCH --error=GWAS/MR_Data/PanUKBB/sex_stratified/logs/format_BS_%j.err

set -e
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript GWAS/MR_Data/PanUKBB/format_panukbb_both_sexes.R
