#!/bin/bash
#SBATCH --job-name=7110
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F7a/logs/%x-%j.out
set -euo pipefail
export MASLD_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
export FIG7_EXEC_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F7a/exec
export FIG7_SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F7a/scripts
export FIG7_OUT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/F7a
export TMPDIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F7a/tmp PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
for f in fig7d_measured_bounds_and_join_structure fig7e_verdict_vocabulary_and_refusal_thresholds fig7f_liver_encoder_against_the_field; do
  "$RS" -e 'set.seed(20260923)' -e "source('$FIG7_SCRIPTS/$f.R')" 2> /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F7a/logs/$f.caption.txt
done
"$RS" -e 'suppressPackageStartupMessages({library(ggplot2);library(patchwork)}); writeLines(capture.output(sessionInfo()), "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F7a/logs/sessionInfo.txt")'
