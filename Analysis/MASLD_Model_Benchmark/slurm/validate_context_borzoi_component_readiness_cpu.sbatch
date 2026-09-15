#!/usr/bin/env bash
#SBATCH --job-name=model-check-981
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONTRACT="${ROOT}/config/artifacts/models/context_borzoi/component_interface_readiness_20260825.json"
TARGET="${ROOT}/executions/model-check-981-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite component-readiness validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"

for pair in \
  "8195a56a52468077d089942b1cd7dec93ce1b728216a94208d07c0c13e231f70 config/artifacts/models/context_borzoi/component_interface_readiness_20260825.json" \
  "b4889e00e4dcd6e19d9e81544ce1aa6bed683017d66fee9b964c8f03c179e2a0 scripts/audit_context_borzoi_component_readiness.py" \
  "53a3b6fb9bf73368799696962745c0033e19299597632a8d110ee3324898aac2 tests/unit/test_audit_context_borzoi_component_readiness.py"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.unit.test_audit_context_borzoi_component_readiness \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}/src:${ROOT}" python \
  scripts/audit_context_borzoi_component_readiness.py \
  --project-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --output "${STAGE}/readiness/receipt.json" \
  > "${STAGE}/validation/audit.stdout.json"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  config/artifacts/models/context_borzoi/component_interface_readiness_20260825.json \
  scripts/audit_context_borzoi_component_readiness.py \
  tests/unit/test_audit_context_borzoi_component_readiness.py \
  slurm/validate_context_borzoi_component_readiness_cpu.sbatch \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "readiness/receipt.json").read_text(encoding="utf-8"))
stderr = (root / "validation/tests.stderr.txt").read_text(encoding="utf-8")
if (
    receipt.get("status") != "pass_interface_contract_build_blocked"
    or receipt.get("verified_frozen_artifacts") != 8
    or receipt.get("component_interface_schema_ready") is not True
    or receipt.get("complementarity_trigger_passed") is not False
    or receipt.get("component_selection_locked") is not False
    or receipt.get("architecture_built") is not False
    or receipt.get("checkpoint_opened_by_audit") is not False
    or receipt.get("biological_data_read") is not False
    or receipt.get("development_outcomes_read") is not False
    or receipt.get("sealed_features_read") is not False
    or receipt.get("sealed_labels_or_outcomes_read") is not False
    or receipt.get("training_authorized") is not False
    or receipt.get("prediction_authorized") is not False
    or receipt.get("gpu_bundle_registration_authorized") is not False
    or receipt.get("gpu_job_submission_authorized") is not False
    or "Ran 4 tests" not in stderr
    or not stderr.rstrip().endswith("OK")
):
    raise SystemExit("component-readiness validation receipt differs")

freeze_tree(
    root,
    {
        "artifact_class": "context_borzoi_component_interface_readiness",
        "model_id": "context_borzoi",
        "artifact_role": "outcome_blind_interface_contract_only",
        "tests": 4,
        "verified_frozen_artifacts": 8,
        "component_interface_schema_ready": True,
        "complementarity_trigger_passed": False,
        "component_selection_locked": False,
        "architecture_built": False,
        "checkpoint_opened": False,
        "biological_data_read": False,
        "development_outcomes_read": False,
        "sealed_features_read": False,
        "sealed_labels_or_outcomes_read": False,
        "training_authorized": False,
        "prediction_authorized": False,
        "gpu_bundle_registered": False,
        "gpu_job_submitted": False,
        "global_frozen_census_modified": False,
        "global_promotion_gate_modified": False,
        "contract_sha256": receipt["contract_sha256"],
        "status": "passed",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Component-readiness validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
