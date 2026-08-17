#!/bin/bash
#SBATCH --job-name=atacaudit
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_audit_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_audit_%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PACKAGE=${ROOT}/Analysis/Multimodal_Program_Projection/scripts/atac_context_v3
CANDIDATE=${ROOT}/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1
LOG_DIR=${CANDIDATE}/audit_logs
MAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba

if [[ -e "${LOG_DIR}" ]]; then
  echo "Refusing to overwrite candidate audit logs: ${LOG_DIR}" >&2
  exit 1
fi
mkdir -p "${LOG_DIR}"
eval "$(${MAMBA} shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u
cd "${ROOT}"

find "${PACKAGE}" -type f -name '*.py' -print0 \
  | xargs -0 python -m py_compile
for script in "${PACKAGE}"/*.R; do
  Rscript -e 'parse(file=commandArgs(TRUE)[1]); cat("R_PARSE_OK\n")' "${script}"
done
bash "${PACKAGE}/tests/run_tests.sh" 2>&1 | tee "${LOG_DIR}/fixture_tests.log"
git diff --check -- \
  Analysis/ATAC/Integration/scripts/35_atac_integration.py \
  Analysis/Multimodal_Program_Projection/scripts/atac_context_v3
python3 scripts/manuscript/validate_resource_scope.py 2>&1 \
  | tee "${LOG_DIR}/resource_scope.log"
python "${PACKAGE}/19_finalize_candidate_audit.py" \
  --scope-log "${LOG_DIR}/resource_scope.log" \
  --tests-log "${LOG_DIR}/fixture_tests.log"
