#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --job-name=pseudobulkde
#SBATCH --output=Analysis/SingleCell/logs/pseudobulkde_%j.out
#SBATCH --error=Analysis/SingleCell/logs/pseudobulkde_%j.err

set -eo pipefail
mkdir -p Analysis/SingleCell/logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript Analysis/SingleCell/scripts/pseudobulk_de_donorcollapse.R
