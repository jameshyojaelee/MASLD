#!/usr/bin/env bash
#SBATCH --job-name=model-incident-001
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
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
SOURCE=${ROOT}/config/artifacts/incidents/gse268273_evaluator_receipt_exposure_20260824.json
SOURCE_SHA256=2b05f192fe072838b09cb0f8dd55ef22a7ca4000d8e6cc39096ea3b4278e4711
TARGET=${ROOT}/executions/gse268273-evaluator-receipt-incident-20260824
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE268273 incident artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SOURCE_SHA256}" "${SOURCE}" \
  | sha256sum --check --strict
cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m unittest -v \
  tests.unit.test_gse268273_evaluator_receipt_incident \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${PYTHON}" scripts/validate_gse268273_evaluator_receipt_incident.py \
  --source "${SOURCE}" \
  --output "${STAGE}/incident" \
  > "${STAGE}/validation/validator_stdout.json" \
  2> "${STAGE}/validation/validator_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  config/artifacts/incidents/gse268273_evaluator_receipt_exposure_20260824.json \
  scripts/validate_gse268273_evaluator_receipt_incident.py \
  tests/unit/test_gse268273_evaluator_receipt_incident.py \
  slurm/validate_gse268273_evaluator_receipt_incident_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
manifest = verify_frozen_tree(root / "incident")
receipt = json.loads((root / "incident/incident_receipt.json").read_text())
if (
    set(receipt)
    != {
        "path_opened",
        "timestamp",
        "turn",
        "outcomes.tsv_not_opened",
        "labels_not_joined",
        "disposition",
    }
    or receipt["outcomes.tsv_not_opened"] is not True
    or receipt["labels_not_joined"] is not True
    or receipt["disposition"] != "evaluator_receipt_exposed"
    or manifest["metadata"]["clean_or_sealed_champion_eligible"] is not False
):
    raise SystemExit("content-free incident firewall differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse268273_evaluator_receipt_incident_validation",
        "clean_or_sealed_champion_eligible": False,
        "independent_unexposed_rederivation_or_audit_required": True,
        "status": "contained_fail_closed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/incident/ARTIFACTS.json"
