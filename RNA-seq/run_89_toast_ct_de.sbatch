#!/bin/bash
#SBATCH --partition=cpu --cpus-per-task=8 --mem=200G --time=48:00:00
#SBATCH --job-name=89_TOAST --output=RNA-seq/logs/89_TOAST_%j.out --error=RNA-seq/logs/89_TOAST_%j.err
set -eo pipefail; mkdir -p RNA-seq/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript RNA-seq/89_toast_celltype_de.R
