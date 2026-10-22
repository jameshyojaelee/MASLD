#!/bin/bash
#SBATCH --job-name=figure7_review
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=logs/figure7_review_%j.out
#SBATCH --error=logs/figure7_review_%j.err
set -euo pipefail
export FIG34_OUT="$PWD/figures/candidates/figure7-native-review-proposal-$(date -u +%Y%m%dT%H%M%SZ)_${SLURM_JOB_ID}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONNOUSERSITE=1
mkdir -p "$FIG34_OUT/source"
printf '%s\n' "$FIG34_OUT" > "logs/figure7_review_${SLURM_JOB_ID}.root"
cp scripts/analysis/bulk_refit_followup/assemble_figure7_review.py scripts/analysis/bulk_refit_followup/run_figure7_review.sbatch "$FIG34_OUT/source/"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python "$FIG34_OUT/source/assemble_figure7_review.py"
