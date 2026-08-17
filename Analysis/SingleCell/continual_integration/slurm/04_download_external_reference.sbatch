#!/usr/bin/env bash
#SBATCH --job-name=masld_hlica_get
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/logs/hlica_download_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/logs/hlica_download_%j.log

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PIPELINE_ROOT="${PROJECT_ROOT}/Analysis/SingleCell/continual_integration"
PYTHON_BIN="${PROJECT_ROOT}/.agent/scanpy_env/bin/python"
RUN_STAMP="$(date -u +%Y-%m-%dT%H%M%SZ)"
RUN_DIR="${PIPELINE_ROOT}/candidates/external-reference-${RUN_STAMP}-${SLURM_JOB_ID}"
DATASET_VERSION="05955131-98d8-4be1-b419-efb690a4e258"
EXPECTED_BYTES="5718629476"
URL="https://datasets.cellxgene.cziscience.com/${DATASET_VERSION}.h5ad"
H5AD="${RUN_DIR}/HLiCA_all_cells_${DATASET_VERSION}.h5ad"

mkdir -p "${RUN_DIR}"
printf '%s\n' "${RUN_DIR}" > "${PIPELINE_ROOT}/logs/hlica_download_${SLURM_JOB_ID}.run_dir"
trap 'status=$?; if [[ ${status} -ne 0 ]]; then printf "%s\n" "${status}" > "${RUN_DIR}/FAILED.exit_code"; fi' EXIT
cd "${PROJECT_ROOT}"
export PYTHONNOUSERSITE=1
export PYTHONPATH="${PIPELINE_ROOT}"

sha256sum \
  "${PIPELINE_ROOT}/reference_policy_v1.json" \
  "${PIPELINE_ROOT}/external_reference_registry_v1.tsv" \
  "${PIPELINE_ROOT}/scripts/inspect_external_reference.py" \
  > "${RUN_DIR}/source_sha256.txt"
printf '%s\n' "${URL}" > "${RUN_DIR}/source_url.txt"

curl --fail --location --retry 5 --retry-all-errors \
  --output "${H5AD}.part" "${URL}"
observed_bytes="$(stat --format='%s' "${H5AD}.part")"
if [[ "${observed_bytes}" != "${EXPECTED_BYTES}" ]]; then
  printf 'external H5AD size mismatch: %s != %s\n' "${observed_bytes}" "${EXPECTED_BYTES}" >&2
  exit 1
fi
mv "${H5AD}.part" "${H5AD}"
sha256sum "${H5AD}" > "${RUN_DIR}/h5ad_sha256.txt"

"${PYTHON_BIN}" "${PIPELINE_ROOT}/scripts/inspect_external_reference.py" \
  --h5ad "${H5AD}" \
  --expected-bytes "${EXPECTED_BYTES}" \
  --expected-dataset-version "${DATASET_VERSION}" \
  --output "${RUN_DIR}/inspection.json" \
  2>&1 | tee "${RUN_DIR}/inspection.log"
"${PYTHON_BIN}" -m pip freeze > "${RUN_DIR}/pip_freeze.txt"
touch "${RUN_DIR}/EXTERNAL_REFERENCE_INSPECTION_COMPLETE"
