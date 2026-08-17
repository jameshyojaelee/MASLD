#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 1 ]]; then
  echo "usage: interactive-gpu-single-runner-v1.sh SPEC.json" >&2
  exit 2
fi
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PIPELINE_ROOT="${PROJECT_ROOT}/Analysis/SingleCell/continual_integration"
PYTHON_BIN="${PROJECT_ROOT}/.agent/scanpy_env/bin/python"
SPEC="$(realpath "$1")"
RUN_STAMP="$(date -u +%Y-%m-%dT%H%M%SZ)"
ALLOC_DIR="${PIPELINE_ROOT}/candidates/execution-gpu-single-${RUN_STAMP}-${SLURM_JOB_ID}"
RECORD="${ALLOC_DIR}/replica"
mkdir -p "${RECORD}"
trap 'status=$?; if [[ ${status} -ne 0 ]]; then printf "%s\n" "${status}" > "${ALLOC_DIR}/FAILED.exit_code"; fi' EXIT
export PYTHONNOUSERSITE=1 PYTHONPATH="${PIPELINE_ROOT}" PYTHONUNBUFFERED=1
export MASLD_RESOURCE_CLASS=large MASLD_PYTHON_BIN="${PYTHON_BIN}"
cd "${PROJECT_ROOT}"
sha256sum "$0" "${SPEC}" "${SPEC}.source-lock.json" > "${ALLOC_DIR}/launch_sha256.txt"
"${PYTHON_BIN}" -m pip freeze > "${ALLOC_DIR}/pip_freeze.txt"
nvidia-smi -q > "${ALLOC_DIR}/nvidia_smi.txt"

set +e
CUDA_VISIBLE_DEVICES=0 OMP_NUM_THREADS=8 "${PYTHON_BIN}" \
  "${PIPELINE_ROOT}/scripts/execute_run_spec.py" \
  --spec "${SPEC}" --pipeline-root "${PIPELINE_ROOT}" --record-dir "${RECORD}" \
  > "${ALLOC_DIR}/replica.log" 2>&1
STATUS=$?
set -e
printf '%s\n' "${STATUS}" > "${ALLOC_DIR}/replica.exit_code"
if [[ ${STATUS} -ne 0 ]]; then
  exit "${STATUS}"
fi
touch "${RECORD}/SLURM_COMPLETE" "${ALLOC_DIR}/GPU_SINGLE_COMPLETE"
