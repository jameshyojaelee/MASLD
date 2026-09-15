#!/usr/bin/env bash
#SBATCH --job-name=model-validate-208
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
TARGET=${ROOT}/executions/model-validate-208-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE268273 source admission validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m compileall -q \
  scripts/admit_gse268273_raw_rsem_source_integrity.py \
  tests/unit/test_admit_gse268273_raw_rsem_source_integrity.py
"${PYTHON}" -m unittest -v \
  tests.unit.test_admit_gse268273_raw_rsem_source_integrity \
  tests.unit.test_gse268273_evaluator_receipt_incident \
  tests.unit.test_gse268273_raw_campaign \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/admit_gse268273_raw_rsem_source_integrity.py \
  tests/unit/test_admit_gse268273_raw_rsem_source_integrity.py \
  config/artifacts/incidents/gse268273_evaluator_receipt_exposure_20260824.json \
  scripts/validate_gse268273_evaluator_receipt_incident.py \
  tests/unit/test_gse268273_evaluator_receipt_incident.py \
  tests/unit/test_gse268273_raw_campaign.py \
  slurm/validate_gse268273_raw_source_integrity_admission_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    root,
    {
        "artifact_class": "gse268273_raw_source_integrity_admission_validation",
        "source_integrity_only": True,
        "clean_or_sealed_champion_eligible": False,
        "independent_unexposed_rederivation_or_audit_required": True,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
