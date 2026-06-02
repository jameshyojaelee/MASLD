#!/bin/bash
#SBATCH --partition=cpu --cpus-per-task=2 --mem=10G --time=48:00:00
#SBATCH --job-name=326_E1
#SBATCH --output=Analysis/SingleCell/logs/326_E1_%j.out
#SBATCH --error=Analysis/SingleCell/logs/326_E1_%j.err
set -eo pipefail
mkdir -p Analysis/SingleCell/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript Analysis/SingleCell/scripts/326_hep_metabolic_macrophage_ccc.R
