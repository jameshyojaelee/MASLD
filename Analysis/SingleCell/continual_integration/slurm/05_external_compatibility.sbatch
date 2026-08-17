#!/usr/bin/env bash
#SBATCH --job-name=masld_hlica_audit
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/logs/hlica_audit_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/logs/hlica_audit_%j.log

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PIPELINE_ROOT="${PROJECT_ROOT}/Analysis/SingleCell/continual_integration"
PYTHON_BIN="${PROJECT_ROOT}/.agent/scanpy_env/bin/python"
RUN_STAMP="$(date -u +%Y-%m-%dT%H%M%SZ)"
RUN_DIR="${PIPELINE_ROOT}/candidates/external-compatibility-${RUN_STAMP}-${SLURM_JOB_ID}"
EXTERNAL_H5AD="${PIPELINE_ROOT}/candidates/external-reference-2026-08-13T234837Z-20212375/HLiCA_all_cells_05955131-98d8-4be1-b419-efb690a4e258.h5ad"
PREPARED_H5AD="${PIPELINE_ROOT}/candidates/prepared-2026-08-13T230932Z-20208265/prepared_counts.h5ad"

mkdir -p "${RUN_DIR}"
printf '%s\n' "${RUN_DIR}" > "${PIPELINE_ROOT}/logs/hlica_audit_${SLURM_JOB_ID}.run_dir"
trap 'status=$?; if [[ ${status} -eq 0 ]]; then touch "${RUN_DIR}/EXTERNAL_COMPATIBILITY_COMPLETE"; else printf "%s\n" "${status}" > "${RUN_DIR}/FAILED.exit_code"; fi' EXIT
cd "${PROJECT_ROOT}"
export PYTHONNOUSERSITE=1
export PYTHONPATH="${PIPELINE_ROOT}"
export PYTHONUNBUFFERED=1
export PYTHONFAULTHANDLER=1

sha256sum \
  "${PIPELINE_ROOT}/config_v1.json" \
  "${PIPELINE_ROOT}/reference_policy_v1.json" \
  "${PIPELINE_ROOT}/external_label_map_v1.tsv" \
  "${PIPELINE_ROOT}/masld_cl/external_compatibility.py" \
  "${PIPELINE_ROOT}/scripts/audit_external_compatibility.py" \
  > "${RUN_DIR}/source_sha256.txt"
sha256sum "${EXTERNAL_H5AD}" "${PREPARED_H5AD}" > "${RUN_DIR}/input_sha256.txt"

"${PYTHON_BIN}" -m py_compile \
  "${PIPELINE_ROOT}/masld_cl/external_compatibility.py" \
  "${PIPELINE_ROOT}/scripts/audit_external_compatibility.py"
"${PYTHON_BIN}" -m unittest discover -v "${PIPELINE_ROOT}/tests" \
  2>&1 | tee "${RUN_DIR}/unit_tests.log"
"${PYTHON_BIN}" "${PIPELINE_ROOT}/scripts/audit_external_compatibility.py" \
  --config "${PIPELINE_ROOT}/config_v1.json" \
  --policy "${PIPELINE_ROOT}/reference_policy_v1.json" \
  --label-map "${PIPELINE_ROOT}/external_label_map_v1.tsv" \
  --external-h5ad "${EXTERNAL_H5AD}" \
  --prepared-h5ad "${PREPARED_H5AD}" \
  --output "${RUN_DIR}/external_compatibility.json" \
  2>&1 | tee "${RUN_DIR}/external_compatibility.log"
"${PYTHON_BIN}" -m pip freeze > "${RUN_DIR}/pip_freeze.txt"
