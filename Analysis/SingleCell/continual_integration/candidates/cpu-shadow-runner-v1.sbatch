#!/usr/bin/env bash
#SBATCH --job-name=masld_cl_cpu_shadow
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=192G
#SBATCH --time=72:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/logs/cpu_shadow_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/logs/cpu_shadow_%j.log

set -euo pipefail
if [[ $# -ne 1 ]]; then
  echo "usage: sbatch cpu-shadow-runner-v1.sbatch EXECUTION_SPEC.json" >&2
  exit 2
fi
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PIPELINE_ROOT="${PROJECT_ROOT}/Analysis/SingleCell/continual_integration"
PYTHON_BIN="${PROJECT_ROOT}/.agent/scanpy_env/bin/python"
SPEC_PATH="$(realpath "$1")"
RUN_STAMP="$(date -u +%Y-%m-%dT%H%M%SZ)"
RECORD_DIR="${PIPELINE_ROOT}/candidates/execution-cpu-shadow-${RUN_STAMP}-${SLURM_JOB_ID}"
mkdir -p "${RECORD_DIR}"
trap 'status=$?; if [[ ${status} -ne 0 ]]; then printf "%s\n" "${status}" > "${RECORD_DIR}/FAILED.exit_code"; fi' EXIT
export PYTHONNOUSERSITE=1 PYTHONPATH="${PIPELINE_ROOT}" PYTHONUNBUFFERED=1
export MASLD_RESOURCE_CLASS=large MASLD_PYTHON_BIN="${PYTHON_BIN}"
cd "${PROJECT_ROOT}"
sha256sum "$0" > "${RECORD_DIR}/runner_sha256.txt"
"${PYTHON_BIN}" -m pip freeze > "${RECORD_DIR}/pip_freeze.txt"
lscpu > "${RECORD_DIR}/lscpu.txt"
free -h > "${RECORD_DIR}/memory.txt"
"${PYTHON_BIN}" "${PIPELINE_ROOT}/scripts/execute_run_spec.py" \
  --spec "${SPEC_PATH}" --pipeline-root "${PIPELINE_ROOT}" --record-dir "${RECORD_DIR}"
touch "${RECORD_DIR}/SLURM_COMPLETE"
