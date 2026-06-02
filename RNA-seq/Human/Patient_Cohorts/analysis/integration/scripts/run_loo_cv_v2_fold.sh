#!/bin/bash
#SBATCH --partition=cpu,io
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --qos=nslab
#SBATCH --job-name=dream-loocv

set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts

echo "=== LOO-CV v2 fold: HELD_OUT=${HELD_OUT} ==="
echo "Started: $(date)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"

Rscript dream_loo_cv_v2.R

echo "Finished: $(date)"
