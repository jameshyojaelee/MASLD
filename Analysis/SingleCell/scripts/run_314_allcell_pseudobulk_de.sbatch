#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --job-name=314_allcell_pb
#SBATCH --output=Analysis/SingleCell/logs/314_allcell_pb_%j.out
#SBATCH --error=Analysis/SingleCell/logs/314_allcell_pb_%j.err

set -eo pipefail
mkdir -p Analysis/SingleCell/logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript Analysis/SingleCell/scripts/314_allcell_pseudobulk_de.R
