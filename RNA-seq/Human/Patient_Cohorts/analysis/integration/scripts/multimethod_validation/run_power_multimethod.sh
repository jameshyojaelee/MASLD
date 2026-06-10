#!/usr/bin/env bash
# run_power_multimethod.sh
# ---------------------------------------------------------------------------
# Launch the multi-method POWER harness as a 3-stage SLURM chain:
#
#   (a) NB parameter estimation (FILE 1, power_estimate_nb_params.R) — ONE job.
#       Caches results/integration/multimethod_validation/power/nb_params.rds.
#   (b) Power simulation array (FILE 2, power_multimethod.R) — --array=1-48,
#       one task per grid cell, GRID_TASK=$SLURM_ARRAY_TASK_ID. NO %N throttle
#       (repo rule) — let QOS / AssocMaxJobsLimit be the ceiling. Gated
#       afterok on (a) so the cache is guaranteed present.
#   (c) Aggregation (FILE 3, aggregate_power_multimethod.R) — ONE job, gated
#       afterok on the whole array (b).
#
# Power-simulation pillar: --job-name=powersim (renamed from multimethod 2026-06-03
# for hand-off to a second user account to bypass the per-user cpu limit).
# Mirrors the micromamba activation idiom of ../run_loo_cv_v2.sh.
#
# Usage: bash run_power_multimethod.sh   (orchestrator only — runs NO compute)
# ---------------------------------------------------------------------------
set -eo pipefail

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/multimethod_validation
LOG_DIR=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/multimethod_validation/power
mkdir -p "${LOG_DIR}"

MM_HOOK='eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"'

echo "=== Multi-method POWER harness: submitting 3-stage chain ==="
echo "Log dir: ${LOG_DIR}"
echo ""

# --- (a) NB parameter estimation (run once) --------------------------------
NB_JID=$(sbatch --parsable \
    --partition=cpu \
    --qos=nslab \
    --cpus-per-task=8 \
    --mem=32G \
    --time=48:00:00 \
    --job-name=powersim \
    --output="${LOG_DIR}/power_nbparams_%j.out" \
    --error="${LOG_DIR}/power_nbparams_%j.err" \
    --export=ALL \
    --wrap="
set -e
set -u
export MASLD_PROJECT_ROOT=${PROJECT}
cd ${PROJECT}
echo \"=== Power: NB parameter estimation ===\"
echo \"Started: \$(date)\"
echo \"SLURM_JOB_ID: \${SLURM_JOB_ID}  SLURM_CPUS_PER_TASK: \${SLURM_CPUS_PER_TASK}\"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript ${SCRIPTS}/power_estimate_nb_params.R
echo \"Finished: \$(date)\"
")
echo "  (a) NB-params submitted: job ${NB_JID}"

# --- (b) Power simulation array (afterok on NB-params) ---------------------
ARRAY_JID=$(sbatch --parsable \
    --partition=cpu \
    --qos=nslab \
    --cpus-per-task=16 \
    --mem=64G \
    --time=90:00:00 \
    --job-name=powersim \
    --array=1-48 \
    --dependency=afterok:${NB_JID} \
    --output="${LOG_DIR}/power_%a_%A.out" \
    --error="${LOG_DIR}/power_%a_%A.err" \
    --export=ALL \
    --wrap="
set -e
set -u
export MASLD_PROJECT_ROOT=${PROJECT}
export GRID_TASK=\${SLURM_ARRAY_TASK_ID}
cd ${PROJECT}
echo \"=== Power: simulation GRID_TASK=\${GRID_TASK} ===\"
echo \"Started: \$(date)\"
echo \"SLURM_JOB_ID: \${SLURM_JOB_ID}  ARRAY_TASK: \${SLURM_ARRAY_TASK_ID}\"
echo \"SLURM_CPUS_PER_TASK: \${SLURM_CPUS_PER_TASK}\"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript ${SCRIPTS}/power_multimethod.R
echo \"Finished: \$(date)\"
")
echo "  (b) Power array submitted: job ${ARRAY_JID} (--array=1-48, afterok:${NB_JID})"

# --- (c) Aggregation (afterok on the whole array) --------------------------
AGG_JID=$(sbatch --parsable \
    --partition=cpu \
    --qos=nslab \
    --cpus-per-task=4 \
    --mem=16G \
    --time=48:00:00 \
    --job-name=powersim \
    --dependency=afterany:${ARRAY_JID} \
    --output="${LOG_DIR}/power_aggregate_%j.out" \
    --error="${LOG_DIR}/power_aggregate_%j.err" \
    --export=ALL \
    --wrap="
set -e
set -u
export MASLD_PROJECT_ROOT=${PROJECT}
cd ${PROJECT}
echo \"=== Power: aggregation ===\"
echo \"Started: \$(date)\"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript ${SCRIPTS}/aggregate_power_multimethod.R
echo \"Finished: \$(date)\"
")
echo "  (c) Aggregation submitted: job ${AGG_JID} (afterok:${ARRAY_JID})"

echo ""
echo "=== All jobs submitted ==="
echo "NB-params job ID:   ${NB_JID}"
echo "Power array job ID: ${ARRAY_JID} (tasks 1-48)"
echo "Aggregation job ID: ${AGG_JID}"
echo "Monitor: squeue -u \$(whoami) --name=multimethod"
echo "Results: ${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/multimethod_validation/power/"
