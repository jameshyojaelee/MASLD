#!/usr/bin/env bash
#SBATCH --job-name=masld_cl_small
#SBATCH --partition=gpu
#SBATCH --qos=nslab
#SBATCH --gres=gpu:b6k:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/logs/small_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/logs/small_%j.log

set -euo pipefail
if [[ $# -ne 1 ]]; then
  echo "usage: sbatch 20_fit_small.sbatch EXECUTION_SPEC.json" >&2
  exit 2
fi
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PIPELINE_ROOT="${PROJECT_ROOT}/Analysis/SingleCell/continual_integration"
PYTHON_BIN="${PROJECT_ROOT}/.agent/scanpy_env/bin/python"
SPEC_PATH="$(realpath "$1")"
RUN_STAMP="$(date -u +%Y-%m-%dT%H%M%SZ)"
RECORD_DIR="${PIPELINE_ROOT}/candidates/execution-${RUN_STAMP}-${SLURM_JOB_ID}"
mkdir -p "${RECORD_DIR}"
trap 'status=$?; if [[ ${status} -ne 0 ]]; then printf "%s\n" "${status}" > "${RECORD_DIR}/FAILED.exit_code"; fi' EXIT
export PYTHONNOUSERSITE=1 PYTHONPATH="${PIPELINE_ROOT}" PYTHONUNBUFFERED=1
export MASLD_RESOURCE_CLASS=small MASLD_PYTHON_BIN="${PYTHON_BIN}"
cd "${PROJECT_ROOT}"
"${PYTHON_BIN}" -m pip freeze > "${RECORD_DIR}/pip_freeze.txt"
nvidia-smi -q > "${RECORD_DIR}/nvidia_smi.txt"
"${PYTHON_BIN}" "${PIPELINE_ROOT}/scripts/execute_run_spec.py" \
  --spec "${SPEC_PATH}" --pipeline-root "${PIPELINE_ROOT}" --record-dir "${RECORD_DIR}"
touch "${RECORD_DIR}/SLURM_COMPLETE"
