#!/bin/bash
#SBATCH --job-name=B5_9cohort
#SBATCH --partition=bigmem
#SBATCH --mem=500G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/logs/B5_9cohort_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/logs/B5_9cohort_%j.err

mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/logs

source ~/.bashrc
micromamba activate rnaseq
set -eo pipefail

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
Rscript scripts/B5_9cohort_sensitivity.R
