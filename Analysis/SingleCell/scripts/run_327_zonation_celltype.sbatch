#!/bin/bash
#SBATCH --partition=cpu --cpus-per-task=2 --mem=10G --time=48:00:00
#SBATCH --job-name=327_E2 --output=Analysis/SingleCell/logs/327_E2_%j.out --error=Analysis/SingleCell/logs/327_E2_%j.err
set -eo pipefail; mkdir -p Analysis/SingleCell/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript Analysis/SingleCell/scripts/327_zonation_disease_celltype.R
