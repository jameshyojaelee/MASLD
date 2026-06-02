#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=48:00:00
#SBATCH --job-name=337_sexCCC
#SBATCH --output=Analysis/SingleCell/logs/337_sexCCC_%j.out
#SBATCH --error=Analysis/SingleCell/logs/337_sexCCC_%j.err
set -eo pipefail
mkdir -p Analysis/SingleCell/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript Analysis/SingleCell/scripts/337_sex_stratified_ccc.R
