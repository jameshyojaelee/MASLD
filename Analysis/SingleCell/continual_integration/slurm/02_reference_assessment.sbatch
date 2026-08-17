#!/usr/bin/env bash
#SBATCH --job-name=masld_ref_assess
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/logs/reference_assessment_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/logs/reference_assessment_%j.log

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PIPELINE_ROOT="${PROJECT_ROOT}/Analysis/SingleCell/continual_integration"
PYTHON_BIN="${PROJECT_ROOT}/.agent/scanpy_env/bin/python"
RUN_STAMP="$(date -u +%Y-%m-%dT%H%M%SZ)"
RUN_DIR="${PIPELINE_ROOT}/candidates/reference-assessment-${RUN_STAMP}-${SLURM_JOB_ID}"

mkdir -p "${RUN_DIR}"
printf '%s\n' "${RUN_DIR}" > "${PIPELINE_ROOT}/logs/reference_assessment_${SLURM_JOB_ID}.run_dir"
trap 'status=$?; if [[ ${status} -ne 0 ]]; then printf "%s\n" "${status}" > "${RUN_DIR}/FAILED.exit_code"; fi' EXIT
cd "${PROJECT_ROOT}"
export PYTHONNOUSERSITE=1
export PYTHONPATH="${PIPELINE_ROOT}"
export PYTHONUNBUFFERED=1
export PYTHONFAULTHANDLER=1

find "${PIPELINE_ROOT}/masld_cl" "${PIPELINE_ROOT}/scripts" "${PIPELINE_ROOT}/tests" "${PIPELINE_ROOT}/slurm" \
  -type f \( -name '*.py' -o -name '*.sbatch' \) -print0 | sort -z | xargs -0 sha256sum \
  > "${RUN_DIR}/source_sha256.txt"
sha256sum \
  "${PIPELINE_ROOT}/config_v1.json" \
  "${PIPELINE_ROOT}/reference_policy_v1.json" \
  "${PIPELINE_ROOT}/external_reference_registry_v1.tsv" \
  "${PIPELINE_ROOT}/README.md" \
  > "${RUN_DIR}/authority_sha256.txt"
sha256sum docs/{README,PAPER,STATUS,RESULTS,ROADMAP}.md \
  >> "${RUN_DIR}/authority_sha256.txt"

"${PYTHON_BIN}" -m py_compile \
  "${PIPELINE_ROOT}/masld_cl/reference_assessment.py" \
  "${PIPELINE_ROOT}/scripts/assess_reference_candidates.py"
bash -n "${PIPELINE_ROOT}/slurm/02_reference_assessment.sbatch"
"${PYTHON_BIN}" -m unittest discover -v "${PIPELINE_ROOT}/tests" \
  2>&1 | tee "${RUN_DIR}/unit_tests.log"
"${PYTHON_BIN}" "${PIPELINE_ROOT}/scripts/assess_reference_candidates.py" \
  --policy "${PIPELINE_ROOT}/reference_policy_v1.json" \
  --registry "${PIPELINE_ROOT}/external_reference_registry_v1.tsv" \
  --repository-root "${PROJECT_ROOT}" \
  --output "${RUN_DIR}/assessment" \
  2>&1 | tee "${RUN_DIR}/assessment.log"
"${PYTHON_BIN}" scripts/manuscript/validate_resource_scope.py \
  2>&1 | tee "${RUN_DIR}/resource_scope_validation.log"
"${PYTHON_BIN}" -m pip freeze > "${RUN_DIR}/pip_freeze.txt"
touch "${RUN_DIR}/REFERENCE_ASSESSMENT_JOB_COMPLETE"
