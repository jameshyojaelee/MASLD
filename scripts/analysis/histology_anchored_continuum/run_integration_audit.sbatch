#!/usr/bin/env bash
#SBATCH --job-name=hac_integration_audit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/integration-audit-%j.out

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
CORE_DIR="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum"
MOLECULAR_DIR="${CORE_DIR}/molecular_layers"
: "${HAC_INTEGRATION_AUDIT_ROOT:?HAC_INTEGRATION_AUDIT_ROOT is unset}"

if [[ ! -d "${HAC_INTEGRATION_AUDIT_ROOT}" ]]; then
  echo "Audit directory does not exist: ${HAC_INTEGRATION_AUDIT_ROOT}" >&2
  exit 1
fi
if [[ -e "${HAC_INTEGRATION_AUDIT_ROOT}/resource_scope.log" ]]; then
  echo "Refusing to overwrite audit outputs: ${HAC_INTEGRATION_AUDIT_ROOT}" >&2
  exit 1
fi

cd "${PROJECT_ROOT}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

python3 scripts/manuscript/validate_resource_scope.py |
  tee "${HAC_INTEGRATION_AUDIT_ROOT}/resource_scope.log"

"${R_BIN}" --vanilla -e '
files <- list.files("scripts/analysis/histology_anchored_continuum", pattern = "[.]R$",
                    recursive = TRUE, full.names = TRUE)
stopifnot(length(files) >= 20L)
for (path in files) parse(file = path, keep.source = FALSE)
cat("ALL_CONTINUUM_R_FILES_PARSE\t", length(files), "\n", sep = "")
' | tee "${HAC_INTEGRATION_AUDIT_ROOT}/r_parse.log"

"${R_BIN}" --vanilla "${CORE_DIR}/tests/run_tests.R" |
  tee "${HAC_INTEGRATION_AUDIT_ROOT}/core_contract_tests.log"
"${R_BIN}" --vanilla "${MOLECULAR_DIR}/tests/run_tests.R" |
  tee "${HAC_INTEGRATION_AUDIT_ROOT}/molecular_contract_tests.log"

python3 "${CORE_DIR}/validate_integration.py" |
  tee "${HAC_INTEGRATION_AUDIT_ROOT}/artifact_audit.log"

sha256sum \
  config/method_exemptions.json \
  scripts/manuscript/validate_resource_scope.py \
  scripts/analysis/histology_anchored_continuum/validate_integration.py \
  figures/main/fig4_singlecell_programs/CANDIDATE_PANEL_INDEX.tsv \
  > "${HAC_INTEGRATION_AUDIT_ROOT}/integration_inputs.sha256"

"${R_BIN}" --vanilla -e \
  'writeLines(capture.output(sessionInfo()), commandArgs(TRUE)[1])' \
  "${HAC_INTEGRATION_AUDIT_ROOT}/sessionInfo.txt"
python3 -m pip freeze > "${HAC_INTEGRATION_AUDIT_ROOT}/pip-freeze.txt"

echo -e "HAC_INTEGRATION_AUDIT\tPASS"
