#!/bin/bash
#SBATCH --array=1-5
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --job-name=multimethod
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/multimethod_validation/loocv/loocv_%a_%A.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/multimethod_validation/loocv/loocv_%a_%A.err
# ===========================================================================
# Multi-method LOO-CV driver (dual-mode #!/bin/bash — runs the body under bash,
# NOT the dash that `sbatch --wrap` uses, so bash arrays / set -o pipefail work).
#
#   1. Submit mode (run directly):  bash run_loocv_multimethod.sh
#        -> sbatch this file as a 5-task array (one held-out cohort per task),
#           capture the array jobid, submit the aggregator afterany.
#        -> ARRAY_QOS env overrides the array QOS (default nslab; set
#           ARRAY_QOS=interactive to bypass QOSMaxCpuPerUserLimit, max 4 concurrent).
#   2. Array mode (SLURM sets SLURM_ARRAY_TASK_ID): run one fold.
# NO %N throttle (repo rule). Aggregator stays on nslab (dependency-gated).
# ===========================================================================
set -eo pipefail

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/multimethod_validation
LOG_DIR=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/multimethod_validation/loocv
RBIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
SELF=${SCRIPTS}/run_loocv_multimethod.sh

# 5 mega cohorts (config/human_datasets.yaml include_in_mega: true).
COHORTS=(GSE126848 GSE130970 GSE135251 GSE162694 GSE213621)

# ---------------------------------------------------------------------------
# SUBMIT MODE
# ---------------------------------------------------------------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  mkdir -p "${LOG_DIR}"
  ARRAY_QOS="${ARRAY_QOS:-nslab}"
  echo "=== Multi-method LOO-CV: submitting ${#COHORTS[@]}-task array (qos=${ARRAY_QOS}) ==="
  echo "Cohorts: ${COHORTS[*]}"

  ARRAY_JID=$(sbatch --parsable --qos=${ARRAY_QOS} "$SELF")
  echo "  Array submitted: job ${ARRAY_JID} (tasks 1-${#COHORTS[@]}, qos=${ARRAY_QOS})"

  AGG_JID=$(sbatch --parsable \
      --partition=cpu --qos=nslab \
      --cpus-per-task=4 --mem=16G --time=48:00:00 \
      --job-name=multimethod \
      --dependency=afterany:${ARRAY_JID} \
      --output="${LOG_DIR}/loocv_aggregate_%j.out" \
      --error="${LOG_DIR}/loocv_aggregate_%j.err" \
      --wrap "cd ${PROJECT}; export MASLD_PROJECT_ROOT=${PROJECT}; ${RBIN} ${SCRIPTS}/aggregate_loocv_multimethod.R")
  echo "  Aggregation submitted: job ${AGG_JID} (afterany:${ARRAY_JID})"
  echo ""
  echo "Monitor: squeue -u \$(whoami) --name=multimethod"
  echo "Results: ${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/multimethod_validation/loocv/"
  exit 0
fi

# ---------------------------------------------------------------------------
# ARRAY MODE: one fold (runs under bash — array indexing is safe here)
# ---------------------------------------------------------------------------
cd "${PROJECT}"
set -u
export MASLD_PROJECT_ROOT="${PROJECT}"
export HELD_OUT="${COHORTS[$((SLURM_ARRAY_TASK_ID - 1))]}"

echo "=== Multi-method LOO-CV fold: HELD_OUT=${HELD_OUT} ==="
echo "Started: $(date)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}  ARRAY_TASK: ${SLURM_ARRAY_TASK_ID}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK}"
${RBIN} "${SCRIPTS}/loocv_multimethod.R"
echo "Finished: $(date)"
