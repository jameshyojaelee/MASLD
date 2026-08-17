#!/bin/bash
#SBATCH --job-name=atacgenimpl
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_genimpl_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_genimpl_%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PACKAGE=${ROOT}/Analysis/Multimodal_Program_Projection/scripts/atac_context_v3
CANDIDATE=${ROOT}/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1
PENDING=${CANDIDATE}/genetics/.implementation_validation.pending.${SLURM_JOB_ID}
FINAL=${CANDIDATE}/genetics/implementation_validation
MAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba

if [[ -e "${PENDING}" || -e "${FINAL}" || -e "${CANDIDATE}/GENETICS_IMPLEMENTATION_READY" ]]; then
  echo "Refusing to overwrite genetics implementation validation" >&2
  exit 1
fi
if compgen -G "${CANDIDATE}/genetics/.implementation_validation.pending.*" >/dev/null; then
  echo "Unresolved prior genetics implementation validation exists" >&2
  exit 1
fi
mkdir -p "${PENDING}/logs"
NON_GENETIC_SHA=$(sha256sum "${CANDIDATE}/NON_GENETIC_READY" | cut -d' ' -f1)
CANDIDATE_AUDIT_SHA=$(sha256sum "${CANDIDATE}/CANDIDATE_AUDIT_READY" | cut -d' ' -f1)
INTEGRATION_SHA=$(sha256sum "${CANDIDATE}/integration/INTEGRATION_READY" | cut -d' ' -f1)
FIG4F_SHA=$(sha256sum "${CANDIDATE}/integration/FIG4F_READY" | cut -d' ' -f1)
RUNNER_SHA=$(sha256sum "${ROOT}/GWAS/finemapping/src/06_susie_coloc.R" | cut -d' ' -f1)

eval "$(${MAMBA} shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u
cd "${ROOT}"

{
  find "${PACKAGE}" -type f -name '*.py' -print0 | xargs -0 python -m py_compile
  echo PY_COMPILE_PASS
  for script in "${PACKAGE}"/*.R; do
    Rscript -e 'parse(file=commandArgs(TRUE)[1])' "${script}" >/dev/null
  done
  echo R_PARSE_PASS
  git diff --check -- Analysis/Multimodal_Program_Projection/scripts/atac_context_v3
  echo DIFF_CHECK_PASS
} 2>&1 | tee "${PENDING}/logs/static_checks.log"
bash "${PACKAGE}/tests/run_tests.sh" 2>&1 | tee "${PENDING}/logs/fixture_tests.log"
Rscript "${PACKAGE}/tests/smoke_coloc_susie_export.R" \
  "${PACKAGE}/genetics_export_helpers.R" 2>&1 \
  | tee "${PENDING}/logs/coloc_export_fixture.log"
python3 scripts/manuscript/validate_resource_scope.py 2>&1 \
  | tee "${PENDING}/logs/resource_scope.log"
python "${PACKAGE}/20_validate_genetics_implementation.py" \
  --pending-root "${PENDING}" \
  --non-genetic-ready-sha256 "${NON_GENETIC_SHA}" \
  --candidate-audit-ready-sha256 "${CANDIDATE_AUDIT_SHA}" \
  --integration-ready-sha256 "${INTEGRATION_SHA}" \
  --fig4f-ready-sha256 "${FIG4F_SHA}" \
  --runner-before-sha256 "${RUNNER_SHA}"
