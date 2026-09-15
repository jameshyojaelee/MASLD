#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-predlock-259
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${EMBEDDINGS_ROOT:?EMBEDDINGS_ROOT is required}"
: "${EMBEDDINGS_SHA256:?EMBEDDINGS_SHA256 is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
TARGET=${ROOT}/executions/gse256398-transcriptformer-prediction-lock-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE256398 prediction lock.\n' >&2
  exit 1
fi
cd "${ROOT}"
export PYTHONPATH=${ROOT}:${ROOT}/src
"${RUNTIME}/env/bin/python" -m unittest -v \
  tests.unit.test_lock_gse256398_transcriptformer_predictions
"${RUNTIME}/env/bin/python" scripts/lock_gse256398_transcriptformer_predictions.py \
  --embeddings-root "${EMBEDDINGS_ROOT}" \
  --embeddings-sha256 "${EMBEDDINGS_SHA256}" \
  --output "${STAGE}"
"${RUNTIME}/env/bin/python" -m pip freeze --all \
  > "${ROOT}/executions/gse256398-transcriptformer-prediction-lock-${SLURM_JOB_ID}.pip-freeze.txt"
sha256sum "${STAGE}/ARTIFACTS.json"
trap - EXIT
mv "${STAGE}" "${TARGET}"
