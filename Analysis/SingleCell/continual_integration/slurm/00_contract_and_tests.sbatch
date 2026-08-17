#!/usr/bin/env bash
#SBATCH --job-name=masld_cl_preflight
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/logs/preflight_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/logs/preflight_%j.log

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PIPELINE_ROOT="${PROJECT_ROOT}/Analysis/SingleCell/continual_integration"
PYTHON_BIN="${PROJECT_ROOT}/.agent/scanpy_env/bin/python"
RUN_STAMP="$(date -u +%Y-%m-%dT%H%M%SZ)"
RUN_DIR="${PIPELINE_ROOT}/candidates/preflight-${RUN_STAMP}-${SLURM_JOB_ID}"

mkdir -p "${RUN_DIR}"
printf '%s\n' "${RUN_DIR}" > "${PIPELINE_ROOT}/logs/preflight_${SLURM_JOB_ID}.run_dir"
trap 'status=$?; if [[ ${status} -eq 0 ]]; then touch "${RUN_DIR}/PREFLIGHT_COMPLETE"; else printf "%s\n" "${status}" > "${RUN_DIR}/FAILED.exit_code"; fi' EXIT
cd "${PROJECT_ROOT}"
export PYTHONNOUSERSITE=1
export PYTHONPATH="${PIPELINE_ROOT}"
export PYTHONUNBUFFERED=1
export PYTHONFAULTHANDLER=1
module load R/4.4.3

find "${PIPELINE_ROOT}/masld_cl" "${PIPELINE_ROOT}/scripts" "${PIPELINE_ROOT}/tests" "${PIPELINE_ROOT}/slurm" \
  -type f \( -name '*.py' -o -name '*.sbatch' \) -print0 | sort -z | xargs -0 sha256sum \
  > "${RUN_DIR}/source_sha256.txt"
sha256sum Analysis/SingleCell/scripts/lib_donor_collapse.R \
  scripts/manuscript/validate_resource_scope.py docs/{README,PAPER,STATUS,RESULTS,ROADMAP}.md \
  >> "${RUN_DIR}/source_sha256.txt"
sha256sum "${PIPELINE_ROOT}/config_v1.json" > "${RUN_DIR}/config_file_sha256.txt"
sha256sum "${PIPELINE_ROOT}/README.md" > "${RUN_DIR}/readme_sha256.txt"

"${PYTHON_BIN}" -m unittest discover -v "${PIPELINE_ROOT}/tests" \
  2>&1 | tee "${RUN_DIR}/unit_tests.log"
bash -n "${PIPELINE_ROOT}"/slurm/*.sbatch

"${PYTHON_BIN}" "${PIPELINE_ROOT}/scripts/build_contract.py" \
  --config "${PIPELINE_ROOT}/config_v1.json" \
  --output "${RUN_DIR}/contract" \
  2>&1 | tee "${RUN_DIR}/contract.log"

"${PYTHON_BIN}" "${PIPELINE_ROOT}/scripts/smoke_scvi_port.py" \
  2>&1 | tee "${RUN_DIR}/scvi_port_smoke.log"

"${PYTHON_BIN}" "${PIPELINE_ROOT}/scripts/audit_frozen_stage_roster.py" \
  --config "${PIPELINE_ROOT}/config_v1.json" \
  --output "${RUN_DIR}/frozen_stage_roster_audit.json" \
  2>&1 | tee "${RUN_DIR}/frozen_stage_roster_audit.log"

Rscript -e 'source("Analysis/SingleCell/scripts/lib_donor_collapse.R"); contract <- data.table::fread(commandArgs(TRUE)[1]); canonical <- contract[dataset == "Liver_Atlas", library_id]; x <- build_srr_to_donor_map(getwd()); liver <- x[names(x) %in% canonical]; stopifnot(length(liver) == 42L, data.table::uniqueN(liver) == 19L)' \
  "${RUN_DIR}/contract/library_manifest.tsv" \
  2>&1 | tee "${RUN_DIR}/liver_donor_map_test.log"

"${PYTHON_BIN}" scripts/manuscript/validate_resource_scope.py \
  2>&1 | tee "${RUN_DIR}/resource_scope_validation.log"

"${PYTHON_BIN}" -m pip freeze > "${RUN_DIR}/pip_freeze.txt"
