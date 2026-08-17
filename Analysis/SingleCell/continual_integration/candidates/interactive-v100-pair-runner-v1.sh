#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo "usage: interactive-v100-pair-runner-v1.sh SPEC_A.json SPEC_B.json" >&2
  exit 2
fi
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PIPELINE_ROOT="${PROJECT_ROOT}/Analysis/SingleCell/continual_integration"
PYTHON_BIN="${PROJECT_ROOT}/.agent/scanpy_env/bin/python"
SPEC_A="$(realpath "$1")"
SPEC_B="$(realpath "$2")"
RUN_STAMP="$(date -u +%Y-%m-%dT%H%M%SZ)"
ALLOC_DIR="${PIPELINE_ROOT}/candidates/execution-v100-pair-${RUN_STAMP}-${SLURM_JOB_ID}"
RECORD_A="${ALLOC_DIR}/replica_a"
RECORD_B="${ALLOC_DIR}/replica_b"
mkdir -p "${RECORD_A}" "${RECORD_B}"
trap 'status=$?; if [[ ${status} -ne 0 ]]; then printf "%s\n" "${status}" > "${ALLOC_DIR}/FAILED.exit_code"; fi' EXIT
export PYTHONNOUSERSITE=1 PYTHONPATH="${PIPELINE_ROOT}" PYTHONUNBUFFERED=1
export MASLD_RESOURCE_CLASS=large MASLD_PYTHON_BIN="${PYTHON_BIN}"
cd "${PROJECT_ROOT}"
sha256sum "$0" "${SPEC_A}" "${SPEC_A}.source-lock.json" \
  "${SPEC_B}" "${SPEC_B}.source-lock.json" > "${ALLOC_DIR}/launch_sha256.txt"
"${PYTHON_BIN}" -m pip freeze > "${ALLOC_DIR}/pip_freeze.txt"
nvidia-smi -q > "${ALLOC_DIR}/nvidia_smi.txt"

CUDA_VISIBLE_DEVICES=0 OMP_NUM_THREADS=8 "${PYTHON_BIN}" \
  "${PIPELINE_ROOT}/scripts/execute_run_spec.py" \
  --spec "${SPEC_A}" --pipeline-root "${PIPELINE_ROOT}" --record-dir "${RECORD_A}" \
  > "${ALLOC_DIR}/replica_a.log" 2>&1 &
PID_A=$!
CUDA_VISIBLE_DEVICES=1 OMP_NUM_THREADS=8 "${PYTHON_BIN}" \
  "${PIPELINE_ROOT}/scripts/execute_run_spec.py" \
  --spec "${SPEC_B}" --pipeline-root "${PIPELINE_ROOT}" --record-dir "${RECORD_B}" \
  > "${ALLOC_DIR}/replica_b.log" 2>&1 &
PID_B=$!

set +e
wait "${PID_A}"
STATUS_A=$?
wait "${PID_B}"
STATUS_B=$?
set -e
printf '%s\n' "${STATUS_A}" > "${ALLOC_DIR}/replica_a.exit_code"
printf '%s\n' "${STATUS_B}" > "${ALLOC_DIR}/replica_b.exit_code"
if [[ ${STATUS_A} -ne 0 || ${STATUS_B} -ne 0 ]]; then
  exit 1
fi
touch "${RECORD_A}/SLURM_COMPLETE" "${RECORD_B}/SLURM_COMPLETE"
touch "${ALLOC_DIR}/V100_PAIR_COMPLETE"
