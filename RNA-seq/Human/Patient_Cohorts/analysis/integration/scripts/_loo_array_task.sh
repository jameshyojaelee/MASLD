#!/usr/bin/env bash
# ------------------------------------------------------------------
# _loo_array_task.sh
# Helper array-task body invoked by run_mash_loo_sweep.sh.
#
# This script is launched as a SLURM array. Each task selects a single
# held-out cohort from the COHORTS_CSV env var (comma-joined) using
# $SLURM_ARRAY_TASK_ID, then runs dream_loo_cv_contrasts.R with the
# (CONTRAST, MASH_DEF, HELD_OUT) triple as environment variables.
#
# Required env vars (exported by the orchestrator via --export):
#   CONTRAST     mash_vs_masl | mash_vs_healthy
#   MASH_DEF     borderline_grouped | strict
#   COHORTS_CSV  comma-joined cohort IDs (index = SLURM_ARRAY_TASK_ID)
# ------------------------------------------------------------------
set -euo pipefail

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
INT=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration
SCRIPTS=${INT}/scripts

# Activate environment
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

# Split COHORTS_CSV into an array and select the held-out cohort
IFS=',' read -r -a cohorts <<<"${COHORTS_CSV}"
if (( SLURM_ARRAY_TASK_ID >= ${#cohorts[@]} )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} out of range (n=${#cohorts[@]})" >&2
    exit 1
fi
HELD_OUT="${cohorts[$SLURM_ARRAY_TASK_ID]}"
export HELD_OUT CONTRAST MASH_DEF

echo "=== LOO array task started: $(date) ==="
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "SLURM_ARRAY_JOB_ID: ${SLURM_ARRAY_JOB_ID}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK:-NA}"
echo "CONTRAST: ${CONTRAST}"
echo "MASH_DEF: ${MASH_DEF}"
echo "HELD_OUT: ${HELD_OUT}"
echo ""

cd "${SCRIPTS}"
Rscript dream_loo_cv_contrasts.R

echo ""
echo "=== LOO array task finished: $(date) ==="
