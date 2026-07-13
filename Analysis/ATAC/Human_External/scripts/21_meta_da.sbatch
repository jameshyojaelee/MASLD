#!/bin/bash
#SBATCH --job-name=limma
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/ATAC/Human_External/logs/metada_%j.out
#SBATCH --error=Analysis/ATAC/Human_External/logs/metada_%j.err
set -eo pipefail
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
Rscript Analysis/ATAC/Human_External/scripts/21_meta_da.R
