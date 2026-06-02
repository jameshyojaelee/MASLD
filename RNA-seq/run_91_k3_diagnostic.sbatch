#!/bin/bash
#SBATCH --partition=cpu --cpus-per-task=4 --mem=32G --time=48:00:00
#SBATCH --job-name=masld_k3diag --output=RNA-seq/logs/masld_k3diag_%j.out --error=RNA-seq/logs/masld_k3diag_%j.err
set -eo pipefail
mkdir -p RNA-seq/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript RNA-seq/91_nmf_k3_diagnostic.R
