#!/bin/bash
#SBATCH --partition=cpu --qos=interactive
#SBATCH --cpus-per-task=8 --mem=120G --time=48:00:00
#SBATCH --job-name=toast_F_trans
#SBATCH --output=RNA-seq/logs/toast_F_trans_%j.out
#SBATCH --error=RNA-seq/logs/toast_F_trans_%j.err
set -eo pipefail
mkdir -p RNA-seq/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

Rscript RNA-seq/scripts/run_toast_celltype_de.R
