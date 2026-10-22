#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
set -euo pipefail
D=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/refit-adoption-review-20261011T151011Z
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript "$D/scripts/review_refit.R" "$1" "$D/out"
