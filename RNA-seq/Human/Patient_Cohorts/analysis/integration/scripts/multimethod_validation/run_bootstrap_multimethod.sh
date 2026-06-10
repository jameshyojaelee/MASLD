#!/bin/bash
#SBATCH --array=1-500
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --job-name=multimethod
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/multimethod_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/multimethod_%A_%a.err
# ===========================================================================
# Multi-method DE validation harness — Pillar A bootstrap driver.
#
# Array of 500 cohort-stratified subsample iterations; each array task runs ALL
# active methods (dream, deseq2; metafor gated) on the SAME resampled split.
# NO %N array throttle — launch full capacity, let QOS / AssocMaxJobsLimit be
# the ceiling (repo rule).
#
# TWO MODES (the SBATCH headers above are only consumed in array mode):
#   1. Submit mode  (run directly, NOT under SLURM):
#        bash run_bootstrap_multimethod.sh
#      -> sbatch --parsable this file as the array, capture the array jobid,
#         then sbatch the aggregator with --dependency=afterok:<arrayjobid>.
#         Echoes both jobids. Does NOT run any compute locally.
#   2. Array mode   (SLURM sets SLURM_ARRAY_TASK_ID):
#        the body runs one iteration via bootstrap_multimethod.R.
# ===========================================================================
set -eo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPT_DIR="RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/multimethod_validation"
LOG_DIR="RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs"
SELF="$PROJECT_ROOT/$SCRIPT_DIR/run_bootstrap_multimethod.sh"

# ---------------------------------------------------------------------------
# SUBMIT MODE: not running inside a SLURM array task -> orchestrate the chain.
# ---------------------------------------------------------------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  cd "$PROJECT_ROOT"
  mkdir -p "$PROJECT_ROOT/$LOG_DIR"

  # 1. Submit the 500-task array (this same file; headers above apply).
  ARRAY_JID=$(sbatch --parsable "$SELF")
  echo "Submitted bootstrap array job: ${ARRAY_JID} (--array=1-500)"

  # 2. Submit the aggregator gated on the array completing OK.
  AGG_JID=$(sbatch --parsable \
    --partition=cpu --qos=nslab \
    --cpus-per-task=4 --mem=32G --time=48:00:00 \
    --job-name=multimethod \
    --output="${LOG_DIR}/multimethod_aggregate_%j.out" \
    --error="${LOG_DIR}/multimethod_aggregate_%j.err" \
    --dependency="afterany:${ARRAY_JID}" \
    --wrap "cd ${PROJECT_ROOT}; \
            export MASLD_PROJECT_ROOT=${PROJECT_ROOT}; \
            /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript ${SCRIPT_DIR}/aggregate_bootstrap_multimethod.R")
  echo "Submitted aggregator job:     ${AGG_JID} (afterok:${ARRAY_JID})"
  exit 0
fi

# ---------------------------------------------------------------------------
# ARRAY MODE: one iteration. Mirror run_bootstrap_10k.sh micromamba idiom.
# ---------------------------------------------------------------------------
cd "$PROJECT_ROOT"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

export MASLD_PROJECT_ROOT="$PROJECT_ROOT"
export ITER=$SLURM_ARRAY_TASK_ID

Rscript "${SCRIPT_DIR}/bootstrap_multimethod.R"
