#!/usr/bin/env bash
#SBATCH --job-name=masld_cl_prepare
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/logs/prepare_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/logs/prepare_%j.log

set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo "usage: sbatch 10_prepare_large.sbatch CONTRACT_LOCK.json LIBRARY_MANIFEST.tsv" >&2
  exit 2
fi

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PIPELINE_ROOT="${PROJECT_ROOT}/Analysis/SingleCell/continual_integration"
PYTHON_BIN="${PROJECT_ROOT}/.agent/scanpy_env/bin/python"
CONTRACT_LOCK="$(realpath "$1")"
LIBRARY_MANIFEST="$(realpath "$2")"
RUN_STAMP="$(date -u +%Y-%m-%dT%H%M%SZ)"
RUN_DIR="${PIPELINE_ROOT}/candidates/prepared-${RUN_STAMP}-${SLURM_JOB_ID}"

mkdir -p "${RUN_DIR}"
printf '%s\n' "${RUN_DIR}" > "${PIPELINE_ROOT}/logs/prepare_${SLURM_JOB_ID}.run_dir"
trap 'status=$?; if [[ ${status} -eq 0 ]]; then touch "${RUN_DIR}/PREPARATION_COMPLETE"; else printf "%s\n" "${status}" > "${RUN_DIR}/FAILED.exit_code"; fi' EXIT
cd "${PROJECT_ROOT}"
export PYTHONNOUSERSITE=1
export PYTHONPATH="${PIPELINE_ROOT}"
export PYTHONUNBUFFERED=1
export PYTHONFAULTHANDLER=1

find "${PIPELINE_ROOT}/masld_cl" "${PIPELINE_ROOT}/scripts" "${PIPELINE_ROOT}/tests" "${PIPELINE_ROOT}/slurm" \
  -type f \( -name '*.py' -o -name '*.sbatch' \) -print0 | sort -z | xargs -0 sha256sum \
  > "${RUN_DIR}/source_sha256.txt"
sha256sum Analysis/SingleCell/scripts/lib_donor_collapse.R \
  scripts/manuscript/validate_resource_scope.py docs/{README,PAPER,STATUS,RESULTS,ROADMAP}.md \
  >> "${RUN_DIR}/source_sha256.txt"
sha256sum "${PIPELINE_ROOT}/config_v1.json" > "${RUN_DIR}/config_file_sha256.txt"
sha256sum "${PIPELINE_ROOT}/README.md" > "${RUN_DIR}/readme_sha256.txt"
sha256sum "${CONTRACT_LOCK}" "${LIBRARY_MANIFEST}" > "${RUN_DIR}/input_sha256.txt"
hostname > "${RUN_DIR}/hostname.txt"
"${PYTHON_BIN}" -m pip freeze > "${RUN_DIR}/pip_freeze.txt"

/usr/bin/time -v "${PYTHON_BIN}" "${PIPELINE_ROOT}/scripts/prepare_model_data.py" \
  --config "${PIPELINE_ROOT}/config_v1.json" \
  --contract-lock "${CONTRACT_LOCK}" \
  --library-manifest "${LIBRARY_MANIFEST}" \
  --output "${RUN_DIR}/prepared_counts.h5ad" \
  2>&1 | tee "${RUN_DIR}/prepare.log"
