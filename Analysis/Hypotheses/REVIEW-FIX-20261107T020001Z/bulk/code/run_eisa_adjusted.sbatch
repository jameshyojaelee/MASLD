#!/bin/bash
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --job-name=eisa-correction
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=3-18:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/bulk/logs/eisa_%j.log
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/bulk/code/refit_eisa_adjusted.R --out /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/bulk/eisa_adjusted_v3 --permutations 200
