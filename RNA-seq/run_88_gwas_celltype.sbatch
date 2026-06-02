#!/bin/bash
#SBATCH --partition=cpu --cpus-per-task=2 --mem=10G --time=48:00:00
#SBATCH --job-name=88_F2
#SBATCH --output=RNA-seq/logs/88_F2_%j.out
#SBATCH --error=RNA-seq/logs/88_F2_%j.err
set -eo pipefail
mkdir -p RNA-seq/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript RNA-seq/88_gwas_variant_celltype_architecture.R
