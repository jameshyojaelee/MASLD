#!/usr/bin/env bash
#SBATCH --job-name=bayesprism_missing
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --array=0-2
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/bayesprism_missing_%a_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/bayesprism_missing_%a_%j.err

# BayesPrism (InstaPrism) backfill for 3 human datasets that completed MuSiC
# but never ran BayesPrism (February batch omission).
#
# Peak memory from prior run: ~315 GB (job 13144841). 500G provides safe buffer.
# Runtime: ~52 min for one dataset; 48h limit covers worst-case.
#
# NOTE: Uses R directly (not Rscript) because the Rscript binary in this env
# has a hard-coded R_HOME from the original install location. R itself works fine.
#
# Datasets: GSE126848 (53 samples), PRJNA512027 (185 samples), GSE167523 (98 samples)
# Submit: sbatch run_bayesprism_missing_human.sbatch

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV_PREFIX=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/music_deconv

REF_HUMAN="${PROJECT_ROOT}/Analysis/Deconvolution/reference/reference_human_sce.rds"
BULK_BASE="${PROJECT_ROOT}/Analysis/Deconvolution/bulk"
RESULTS_BASE="${PROJECT_ROOT}/Analysis/Deconvolution/results"
SCRIPTS_DIR="${PROJECT_ROOT}/Analysis/Deconvolution/scripts"

# Use R directly — Rscript has stale hard-coded R_HOME from original install path
R_BIN="${ENV_PREFIX}/bin/R"
export R_LIBS_USER=""
export R_LIBS_SITE="${ENV_PREFIX}/lib/R/library"
export LD_LIBRARY_PATH="${ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

DATASETS=(
  "GSE126848"
  "PRJNA512027"
  "GSE167523"
)

DATASET="${DATASETS[$SLURM_ARRAY_TASK_ID]}"

echo "=============================="
echo "BayesPrism backfill: ${DATASET}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "R binary: ${R_BIN}"
echo "=============================="

COUNTS_FILE="${BULK_BASE}/${DATASET}/${DATASET}_counts.tsv"
RESULT_DIR="${RESULTS_BASE}/${DATASET}"

if [[ ! -f "${COUNTS_FILE}" ]]; then
  echo "ERROR: Missing counts file: ${COUNTS_FILE}" >&2
  exit 1
fi

if [[ ! -f "${REF_HUMAN}" ]]; then
  echo "ERROR: Missing human SCE reference: ${REF_HUMAN}" >&2
  exit 1
fi

# Idempotency guard: skip if output already exists
PROP_OUT="${RESULT_DIR}/${DATASET}_bayesprism_proportions.tsv"
if [[ -f "${PROP_OUT}" ]]; then
  echo "WARNING: Output already exists: ${PROP_OUT}"
  echo "Remove it manually to re-run. Exiting."
  exit 0
fi

mkdir -p "${RESULT_DIR}"

echo "Running 12_run_bayesprism.R..."
"${R_BIN}" --no-save --no-restore --slave \
  --file="${SCRIPTS_DIR}/12_run_bayesprism.R" \
  --args \
  "${REF_HUMAN}" \
  "${COUNTS_FILE}" \
  "${RESULT_DIR}" \
  "${DATASET}" \
  "celltype"

echo "Finished: $(date)"
echo "Outputs:"
ls -lh "${RESULT_DIR}/${DATASET}_bayesprism"* 2>/dev/null || echo "  (no output files found -- check for errors)"
