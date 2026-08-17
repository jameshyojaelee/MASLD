#!/usr/bin/env bash
#SBATCH --job-name=masld_cl_harmony
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/candidates/harmony-2026-08-13T230657Z/slurm.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/candidates/harmony-2026-08-13T230657Z/slurm.log

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PIPELINE_ROOT="${PROJECT_ROOT}/Analysis/SingleCell/continual_integration"
PYTHON_BIN="${PROJECT_ROOT}/.agent/scanpy_env/bin/python"
RUN_ROOT="${PIPELINE_ROOT}/candidates/harmony-2026-08-13T230657Z"
CONTRACT_LOCK="${PIPELINE_ROOT}/candidates/preflight-2026-08-13T224413Z-20206577/contract/contract_lock.json"
trap 'status=$?; if [[ ${status} -eq 0 ]]; then touch "${RUN_ROOT}/HARMONY_COMPLETE"; else printf "%s\n" "${status}" > "${RUN_ROOT}/FAILED.exit_code"; fi' EXIT
cd "${PROJECT_ROOT}"
export PYTHONNOUSERSITE=1 PYTHONPATH="${PIPELINE_ROOT}" PYTHONUNBUFFERED=1
find "${PIPELINE_ROOT}/masld_cl" "${PIPELINE_ROOT}/scripts" "${PIPELINE_ROOT}/tests" "${PIPELINE_ROOT}/slurm" \
  -type f \( -name '*.py' -o -name '*.sbatch' \) -print0 | sort -z | xargs -0 sha256sum \
  > "${RUN_ROOT}/source_sha256.txt"
sha256sum "${PIPELINE_ROOT}/config_v1.json" "${CONTRACT_LOCK}" \
  > "${RUN_ROOT}/input_sha256.txt"
"${PYTHON_BIN}" -m pip freeze > "${RUN_ROOT}/pip_freeze.txt"
hostname > "${RUN_ROOT}/hostname.txt"
"${PYTHON_BIN}" "${PIPELINE_ROOT}/scripts/export_harmony.py" \
  --config "${PIPELINE_ROOT}/config_v1.json" \
  --contract-lock "${CONTRACT_LOCK}" \
  --output "${RUN_ROOT}/bundle" \
  > "${RUN_ROOT}/export.json"
