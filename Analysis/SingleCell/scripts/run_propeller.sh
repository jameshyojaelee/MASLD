#!/bin/bash
#SBATCH --job-name=S4_propeller
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/S4_propeller_%j.out
#SBATCH --error=logs/S4_propeller_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
mkdir -p logs
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

Rscript Analysis/SingleCell/scripts/349_propeller_proportion.R
