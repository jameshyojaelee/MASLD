#!/bin/bash
#SBATCH --job-name=s_v3_M07b
#SBATCH --partition=io
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --array=1-5%5
#SBATCH --output=logs/M07b_%A_%a.out
#SBATCH --error=logs/M07b_%A_%a.err
# ---------------------------------------------------------------------------
# Module 07b -- LOCO posterior-predictive array (5 parallel SLURM tasks).
# Each task holds out ONE cohort (Bril, Chen, Govaere, Hoang, Suppli --
# alphabetic, controlled by LOCO_FOLD=$SLURM_ARRAY_TASK_ID) and refits
# dream M2 + mash on the remaining 4 cohorts.
#
# Wraps: scripts/sex_v3/07b_loco_fold.R
# Env:   LOCO_FOLD = $SLURM_ARRAY_TASK_ID (1..5)
#
# Compute notes:
#   - Single LOCO fold: dream M2 on ~4-cohort subset (~700 samples,
#     ~17k genes). Wall ~40-60 min with 16 CPUs.
#   - %5 throttle keeps all 5 folds running concurrently.
# ---------------------------------------------------------------------------
set -e

SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate mashr_env
set -u

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export VECLIB_MAXIMUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

export LOCO_FOLD=${SLURM_ARRAY_TASK_ID}

echo "=== Module 07b (LOCO_FOLD=${LOCO_FOLD}) start: $(date) ==="
echo "Host: $(hostname)  SLURM_JOB_ID: ${SLURM_JOB_ID:-NA}  ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID:-NA}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}  Mem: 96G  Walltime: 8h"
echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
echo

Rscript 07b_loco_fold.R

echo
echo "=== Module 07b (LOCO_FOLD=${LOCO_FOLD}) end: $(date) ==="
