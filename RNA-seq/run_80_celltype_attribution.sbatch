#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=48:00:00
#SBATCH --job-name=80_ct_attr
#SBATCH --output=RNA-seq/logs/80_celltype_attr_%j.out
#SBATCH --error=RNA-seq/logs/80_celltype_attr_%j.err

set -eo pipefail
mkdir -p RNA-seq/logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript RNA-seq/80_celltype_intrinsic_attribution.R
