#!/bin/bash
#SBATCH --job-name=s_v3_M07c
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/M07c_%j.out
#SBATCH --error=logs/M07c_%j.err
# ---------------------------------------------------------------------------
# Module 07c -- Aggregate Test 1 (07a) + 5 LOCO folds (07b) into final
# calibration_qc.pdf + posterior_predictive_per_cohort.csv. Light job:
# table assembly + ggplot rendering.
#
# Wraps: scripts/sex_v3/07c_aggregate.R
# Depends on: afterok 07a + afterok the 07b array.
# ---------------------------------------------------------------------------
set -e

SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate mashr_env
set -u

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

echo "=== Module 07c (aggregate) start: $(date) ==="
echo "Host: $(hostname)  SLURM_JOB_ID: ${SLURM_JOB_ID:-NA}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}  Mem: 32G  Walltime: 2h"
echo

Rscript 07c_aggregate.R

echo
echo "=== Module 07c end: $(date) ==="
