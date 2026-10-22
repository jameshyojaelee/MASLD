#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=5
#SBATCH --mem=48G
set -euo pipefail
D=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/refit-adoption-review-20261011T151011Z
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript "$D/scripts/review_refit_v2.R" loco "$D/out"
