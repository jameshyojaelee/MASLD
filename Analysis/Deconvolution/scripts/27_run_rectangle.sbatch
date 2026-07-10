#!/usr/bin/env bash
#SBATCH --job-name=rectangle
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=0-4
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/rectangle_%a_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/rectangle_%a_%j.err

# ---------------------------------------------------------------------------
# Rectangle deconvolution (rectanglepy 1.5.0) over the 5 pooled-analysis cohorts.
#
# For each cohort, runs script 26 TWICE:
#   1. coarse reference (obs cell_type)      -> 16-type re-normalized proportions
#   2. fine reference  (obs cell_type_fine)  -> all fine labels + Unknown
#
# Prerequisites (produced by sibling steps; assumed present):
#   - reference/rectangle_ref_coarse.h5ad  (obs cell_type, .X raw counts, var=symbols)
#   - reference/rectangle_ref_fine.h5ad    (obs cell_type_fine)
#   - bulk/{COHORT}/{COHORT}_tpm_symbol.tsv (rows=SYMBOL, cols=samples, TPM)
#
# Submit: sbatch 27_run_rectangle.sbatch
#
# Partition fallbacks when cpu is busy (bypass the nslab 7TB cap, max 4 concurrent):
#   io   + --qos=interactive  (fast lane, mem <= ~96-128G) -- ample for 64G here
#   bigmem + --qos=interactive (only if a run needs >210G; enforces a 500G floor)
# ---------------------------------------------------------------------------

set -eo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
RECT_ENV="${PROJECT_ROOT}/RNA-seq/.mamba/rectangle"

eval "$(micromamba shell hook -s bash)"
micromamba activate "${RECT_ENV}"

set -u

# Deterministic hashing for reproducibility (must be set before python starts).
export PYTHONHASHSEED=0

REF_DIR="${PROJECT_ROOT}/Analysis/Deconvolution/reference"
BULK_BASE="${PROJECT_ROOT}/Analysis/Deconvolution/bulk"
RESULTS_BASE="${PROJECT_ROOT}/Analysis/Deconvolution/results"
SCRIPTS_DIR="${PROJECT_ROOT}/Analysis/Deconvolution/scripts"

REF_COARSE="${REF_DIR}/rectangle_ref_coarse.h5ad"
REF_FINE="${REF_DIR}/rectangle_ref_fine.h5ad"
RUNNER="${SCRIPTS_DIR}/26_run_rectangle.py"

# 5 pooled-analysis cohorts (indexed by SLURM_ARRAY_TASK_ID; no %N throttle).
DATASETS=(
  "GSE126848"
  "GSE135251"
  "GSE130970"
  "GSE213621"
  "GSE162694"
)

DATASET="${DATASETS[$SLURM_ARRAY_TASK_ID]}"
BULK_TPM="${BULK_BASE}/${DATASET}/${DATASET}_tpm_symbol.tsv"
RESULT_DIR="${RESULTS_BASE}/${DATASET}"

echo "=============================="
echo "Rectangle deconvolution: ${DATASET}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Host: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-?}"
echo "Started: $(date)"
echo "=============================="

for f in "${REF_COARSE}" "${REF_FINE}" "${BULK_TPM}"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: missing required input: ${f}" >&2
    exit 1
  fi
done

mkdir -p "${RESULT_DIR}"

echo ""
echo "--- [1/2] COARSE reference (cell_type) ---"
srun --ntasks=1 python "${RUNNER}" \
  "${REF_COARSE}" \
  "${BULK_TPM}" \
  "${RESULT_DIR}" \
  "${DATASET}" \
  cell_type

echo ""
echo "--- [2/2] FINE reference (cell_type_fine) ---"
srun --ntasks=1 python "${RUNNER}" \
  "${REF_FINE}" \
  "${BULK_TPM}" \
  "${RESULT_DIR}" \
  "${DATASET}" \
  cell_type_fine

echo ""
echo "Finished: $(date)"
echo "Outputs:"
ls -lh "${RESULT_DIR}/${DATASET}_rectangle"* 2>/dev/null || echo "  (no output files)"
