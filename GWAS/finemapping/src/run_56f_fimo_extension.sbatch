#!/bin/bash
#SBATCH --job-name=fimo_extension
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/56f_fimo_%j.out
#SBATCH --error=logs/56f_fimo_%j.err
set -e
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate motifbreakr
module load meme/5.5.9
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
mkdir -p logs results/gwas_atac
Rscript src/56f_fimo_extension.R
