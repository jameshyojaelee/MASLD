#!/usr/bin/env bash
#SBATCH --job-name=cbx_sigmatrix
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --time=72:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/cbx_sigmatrix_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/cbx_sigmatrix_%j.err

# CIBERSORTx Step A: Build signature matrices from single-cell references.
# Builds one signature matrix for human and one for mouse.
# These matrices are used as input for the deconvolution step (run_cibersortx_deconv.sbatch).
#
# Resources: 400G bigmem for 21GB human sc_matrix (306k cells x 35k genes).
# Runtime up to 72h (building sig matrix from large sc ref is slow).
#
# Prerequisites:
#   1. Fill in credentials: Analysis/Deconvolution/cibersortx_input/.cibersortx_creds
#   2. Pull Singularity image (one-time, on login node or compute node):
#      module load singularity/4.2.2
#      singularity pull --disable-cache \
#        docker://cibersortx/fractions \
#        Analysis/Deconvolution/cibersortx_input/cibersortx_fractions.sif
#
# Submit: sbatch run_cibersortx_build_sig.sbatch

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
CREDS_FILE="${PROJECT_ROOT}/Analysis/Deconvolution/cibersortx_input/.cibersortx_creds"
SIF="${PROJECT_ROOT}/Analysis/Deconvolution/cibersortx_input/cibersortx_fractions.sif"
INPUT_DIR="${PROJECT_ROOT}/Analysis/Deconvolution/cibersortx_input"
SIG_DIR="${PROJECT_ROOT}/Analysis/Deconvolution/results/cibersortx_sigmatrix"

module load singularity/4.2.2

# Load credentials
if [[ ! -f "${CREDS_FILE}" ]]; then
  echo "ERROR: Credentials file not found: ${CREDS_FILE}" >&2
  echo "  Create it with: CIBERSORTX_USER=email CIBERSORTX_TOKEN=token" >&2
  exit 1
fi
source "${CREDS_FILE}"

if [[ "${CIBERSORTX_USER}" == "your_email@domain.com" ]]; then
  echo "ERROR: Edit ${CREDS_FILE} with your actual credentials before submitting" >&2
  exit 1
fi

if [[ ! -f "${SIF}" ]]; then
  echo "ERROR: Singularity image not found: ${SIF}" >&2
  echo "  Pull it with:" >&2
  echo "    module load singularity/4.2.2" >&2
  echo "    singularity pull --disable-cache docker://cibersortx/fractions ${SIF}" >&2
  exit 1
fi

mkdir -p "${SIG_DIR}/human" "${SIG_DIR}/mouse"

echo "=============================="
echo "CIBERSORTx: Build Signature Matrices"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "=============================="

# --- Human signature matrix ---
echo ""
echo "--- Building HUMAN signature matrix ---"
echo "Input: ${INPUT_DIR}/human_sc_matrix.txt (21GB, 306k cells)"

singularity exec \
  --bind "${INPUT_DIR}:/input:ro" \
  --bind "${SIG_DIR}/human:/output" \
  "${SIF}" CIBERSORTxFractions \
  --username "${CIBERSORTX_USER}" \
  --token "${CIBERSORTX_TOKEN}" \
  --single_cell TRUE \
  --refsample /input/human_sc_matrix.txt \
  --phenoclasses /input/human_phenoclasses.txt \
  --outdir /output \
  --verbose TRUE \
  --mode S

echo "Human sig matrix done: $(date)"
ls -lh "${SIG_DIR}/human/"

# --- Mouse signature matrix ---
echo ""
echo "--- Building MOUSE signature matrix ---"
echo "Input: ${INPUT_DIR}/mouse_sc_matrix.txt (29GB)"

singularity exec \
  --bind "${INPUT_DIR}:/input:ro" \
  --bind "${SIG_DIR}/mouse:/output" \
  "${SIF}" CIBERSORTxFractions \
  --username "${CIBERSORTX_USER}" \
  --token "${CIBERSORTX_TOKEN}" \
  --single_cell TRUE \
  --refsample /input/mouse_sc_matrix.txt \
  --phenoclasses /input/mouse_phenoclasses.txt \
  --outdir /output \
  --verbose TRUE \
  --mode S

echo "Mouse sig matrix done: $(date)"
ls -lh "${SIG_DIR}/mouse/"

echo ""
echo "=============================="
echo "Signature matrix build complete: $(date)"
echo "Human: ${SIG_DIR}/human/CIBERSORTx_SigMatrix.txt"
echo "Mouse: ${SIG_DIR}/mouse/CIBERSORTx_SigMatrix.txt"
echo "Next step: sbatch run_cibersortx_deconv.sbatch"
echo "=============================="
