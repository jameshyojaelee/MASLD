#!/usr/bin/env bash
#SBATCH --job-name=cbx_deconv
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --time=48:00:00
#SBATCH --array=0-7
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/cbx_deconv_%a_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/cbx_deconv_%a_%j.err

# CIBERSORTx Step B: Deconvolve bulk RNA-seq for all 8 datasets.
# Uses prebuilt signature matrices from run_cibersortx_build_sig.sbatch.
#
# Resources: 300G bigmem / 48h per task (loading sig matrix + bulk matrix).
#
# DEPENDENCY: run_cibersortx_build_sig.sbatch must complete first.
# Submit with dependency:
#   SIG_JOB=$(sbatch --parsable run_cibersortx_build_sig.sbatch)
#   sbatch --dependency=afterok:${SIG_JOB} run_cibersortx_deconv.sbatch
# Or after sig matrices are ready:
#   sbatch run_cibersortx_deconv.sbatch

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
CREDS_FILE="${PROJECT_ROOT}/Analysis/Deconvolution/cibersortx_input/.cibersortx_creds"
SIF="${PROJECT_ROOT}/Analysis/Deconvolution/cibersortx_input/cibersortx_fractions.sif"
BULK_BASE="${PROJECT_ROOT}/Analysis/Deconvolution/bulk"
RESULTS_BASE="${PROJECT_ROOT}/Analysis/Deconvolution/results"
SIG_DIR="${PROJECT_ROOT}/Analysis/Deconvolution/results/cibersortx_sigmatrix"

module load singularity/4.2.2

# Load credentials
source "${CREDS_FILE}"

# Dataset array: name|species
DATASETS=(
  "GSE130970|human"
  "GSE135251|human"
  "GSE126848|human"
  "PRJNA512027|human"
  "GSE167523|human"
  "GSE156918|mouse"
  "GSE205974|mouse"
  "inhouse_MCD|mouse"
)

IFS='|' read -r DATASET SPECIES <<< "${DATASETS[$SLURM_ARRAY_TASK_ID]}"

echo "=============================="
echo "CIBERSORTx deconvolution: ${DATASET} (${SPECIES})"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "=============================="

BULK_FILE="${BULK_BASE}/${DATASET}/${DATASET}_counts.tsv"
SIG_MATRIX="${SIG_DIR}/${SPECIES}/CIBERSORTx_SigMatrix.txt"
OUT_DIR="${RESULTS_BASE}/${DATASET}/cibersortx"

if [[ ! -f "${BULK_FILE}" ]]; then
  echo "ERROR: Missing bulk counts: ${BULK_FILE}" >&2
  exit 1
fi

if [[ ! -f "${SIG_MATRIX}" ]]; then
  echo "ERROR: Signature matrix not found: ${SIG_MATRIX}" >&2
  echo "  Run run_cibersortx_build_sig.sbatch first." >&2
  exit 1
fi

# Idempotency guard
OUT_FILE="${OUT_DIR}/CIBERSORTx_Results.txt"
if [[ -f "${OUT_FILE}" ]]; then
  echo "WARNING: Output already exists: ${OUT_FILE}. Skipping."
  exit 0
fi

mkdir -p "${OUT_DIR}"

# CIBERSORTx requires write access to the directory containing the sig matrix
# Use a per-job temp dir for output to avoid collisions in array runs
TMP_SIG_DIR=$(mktemp -d)
cp "${SIG_MATRIX}" "${TMP_SIG_DIR}/"
SIG_BASENAME=$(basename "${SIG_MATRIX}")

echo "Running CIBERSORTx fractions deconvolution..."
singularity exec \
  --bind "${BULK_BASE}/${DATASET}:/bulk:ro" \
  --bind "${TMP_SIG_DIR}:/sigmatrix" \
  --bind "${OUT_DIR}:/output" \
  "${SIF}" CIBERSORTxFractions \
  --username "${CIBERSORTX_USER}" \
  --token "${CIBERSORTX_TOKEN}" \
  --sigmatrix "/sigmatrix/${SIG_BASENAME}" \
  --mixture "/bulk/${DATASET}_counts.tsv" \
  --outdir /output \
  --verbose TRUE \
  --absolute FALSE \
  --perm 100

rm -rf "${TMP_SIG_DIR}"

echo "Finished: $(date)"
echo "Results in ${OUT_DIR}:"
ls -lh "${OUT_DIR}/" 2>/dev/null || echo "  (no output files found)"
