#!/bin/bash
#SBATCH --job-name=liver_lincs_rerun
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/lincs_rerun_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/lincs_rerun_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "=== Running Script 31 (LINCS annotation, bug fix) ==="
echo "Started at $(date)"

Rscript RNA-seq/31_lincs_annotation.R

echo "=== Done at $(date) ==="
