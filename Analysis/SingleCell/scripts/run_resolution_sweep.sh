#!/bin/bash
# S3: Submit chained Leiden resolution sweep -> cross-cohort Progressor test -> aggregate report.
# Outputs job IDs to stdout; status JSON consumed by agent_status/S3_STATUS.json.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
LOG_DIR="${SCRIPT_DIR}/logs/resolution_sweep"
mkdir -p "${LOG_DIR}"

echo "=== S3 resolution sweep chain ==="
echo "Log dir: ${LOG_DIR}"

# ── 360: Leiden sweep (GPU) ─────────────────────────────────────────────
JOB_360=$(sbatch --parsable \
  --job-name=s3_360_sweep \
  --partition=gpu \
  --gres=gpu:1 \
  --cpus-per-task=16 \
  --mem=64G \
  --time=48:00:00 \
  --output="${LOG_DIR}/360_sweep_%j.out" \
  --error="${LOG_DIR}/360_sweep_%j.err" \
  --wrap="export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/lib/python3.12/site-packages/nvidia/cu13/lib && micromamba run -n rapids_singlecell python ${SCRIPT_DIR}/360_resolution_sweep.py")
echo "360 submitted: ${JOB_360}"

# ── 361: Cross-cohort Progressor test (CPU) ─────────────────────────────
JOB_361=$(sbatch --parsable \
  --dependency=afterok:${JOB_360} \
  --job-name=s3_361_xcohort \
  --partition=cpu \
  --cpus-per-task=16 \
  --mem=32G \
  --time=48:00:00 \
  --output="${LOG_DIR}/361_xcohort_%j.out" \
  --error="${LOG_DIR}/361_xcohort_%j.err" \
  --wrap="micromamba run -n rapids_singlecell python ${SCRIPT_DIR}/361_cross_cohort_progressor.py")
echo "361 submitted: ${JOB_361} (depends on ${JOB_360})"

# ── 362: Aggregate report (CPU) ─────────────────────────────────────────
JOB_362=$(sbatch --parsable \
  --dependency=afterok:${JOB_361} \
  --job-name=s3_362_report \
  --partition=cpu \
  --cpus-per-task=4 \
  --mem=16G \
  --time=48:00:00 \
  --output="${LOG_DIR}/362_report_%j.out" \
  --error="${LOG_DIR}/362_report_%j.err" \
  --wrap="micromamba run -n rnaseq Rscript ${SCRIPT_DIR}/362_aggregate_resolution_report.R")
echo "362 submitted: ${JOB_362} (depends on ${JOB_361})"

echo ""
echo "=== Chain submitted ==="
echo "360 (sweep):   ${JOB_360}"
echo "361 (xcohort): ${JOB_361}"
echo "362 (report):  ${JOB_362}"

# Emit IDs for status JSON capture
echo "JOB_IDS:${JOB_360},${JOB_361},${JOB_362}"
