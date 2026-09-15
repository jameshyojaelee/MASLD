#!/usr/bin/env bash
#SBATCH --job-name=model-encoder-field-7f-score
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Resource ask re-derived: 21 models x (4 global + 7 per-study) bootstrap
# distributions of 10,000 float64 draws is ~19 MB; the 10,000 x 102 int16
# multiplicity table is 2 MB; the largest transient is the canonical JSON used
# for the distribution digest (~50 MB of text).  32G is far above need and
# matches the sibling scorers; wall clock is minutes, 4 h requested.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

: "${PRODUCER:?frozen encoder-field producer bundle is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SOURCE_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
SPLIT=${ROOT}/executions/model-data-052-21077595/split
SPLIT_SHA256=10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680
REGISTRY=${ROOT}/config/encoder_field_common_head_50000_registry.json
TARGET=${ROOT}/executions/encoder-field-score-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  status=$?
  if [[ ${status} -ne 0 && -d ${STAGE} && ! -e ${FAILED} ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
  exit "${status}"
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite the encoder-field evaluation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
for pair in \
  "${RUNTIME_SHA256} ${RUNTIME}/ARTIFACTS.json" \
  "${SOURCE_SHA256} ${SOURCE}/ARTIFACTS.json" \
  "${SPLIT_SHA256} ${SPLIT}/ARTIFACTS.json" ; do
  printf '%s  %s\n' ${pair} | sha256sum --check --strict
done
if [[ ! -f ${PRODUCER}/ARTIFACTS.json ]]; then
  printf 'Producer bundle is not frozen.\n' >&2
  exit 1
fi
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
cp "${PRODUCER}/validation/guard_report.json" "${STAGE}/validation/producer_guard_report.json"

"${PYTHON}" scripts/score_encoder_field_common_head_50000.py \
  --source "${SOURCE}" \
  --split "${SPLIT}" \
  --registry "${REGISTRY}" \
  --shards "${PRODUCER}/shards" \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/validation/scorer.stdout.json" \
  2> "${STAGE}/validation/scorer.stderr.txt"

sha256sum \
  scripts/score_encoder_field_common_head_50000.py \
  scripts/fit_predict_encoder_field_common_head_50000.py \
  "${REGISTRY}" \
  "${PRODUCER}/ARTIFACTS.json" \
  "${SOURCE}/ARTIFACTS.json" \
  "${SPLIT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

mv "${STAGE}" "${TARGET}"
printf 'encoder field evaluation: %s\n' "${TARGET}"
