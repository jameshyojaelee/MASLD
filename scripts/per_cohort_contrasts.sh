#!/bin/bash
#SBATCH --job-name=per_cohort_contrasts
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/per_cohort_contrasts_%j.out
#SBATCH --error=logs/per_cohort_contrasts_%j.out

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs

micromamba run -n rnaseq Rscript scripts/per_cohort_contrasts.R
