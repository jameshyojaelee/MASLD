#!/usr/bin/env bash
#SBATCH --job-name=model-check-984
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
CONTRACT="${ROOT}/config/artifacts/models/scbasset/production_readiness_20260825.json"
TARGET="${ROOT}/executions/model-check-984-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite scBasset readiness validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"

for pair in \
  "da4cebe5171414bf10474775beda6a58d721a73cb6b13121a3f6daaa1523d577 config/artifacts/models/scbasset/production_readiness_20260825.json" \
  "fa2f6e9fc3071e4d735faa17b2dc24bcd8beb68cc701c122a7afc7e277f72b41 scripts/audit_scbasset_production_readiness.py" \
  "e1ba6eaa605c41fc2e8a54942261c6558c8fe7b52f86d78d12d24cba7e432519 tests/unit/test_audit_scbasset_production_readiness.py"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.unit.test_audit_scbasset_production_readiness \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}/src:${ROOT}" python \
  scripts/audit_scbasset_production_readiness.py \
  --project-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --output "${STAGE}/readiness/receipt.json" \
  > "${STAGE}/validation/audit.stdout.txt"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  config/artifacts/models/scbasset/production_readiness_20260825.json \
  scripts/audit_scbasset_production_readiness.py \
  tests/unit/test_audit_scbasset_production_readiness.py \
  slurm/validate_scbasset_production_readiness_cpu.sbatch \
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
    receipt.get("status") != "pass_reconciled_development_only_no_gpu_unit"
    or receipt.get("native_observed_atac_role")
    != "outer_training_supervision_only"
    or receipt.get("observed_query_atac_consumed") is not False
    or receipt.get("inductive_observed_multiome_model") is not False
    or receipt.get("five_outer_fold_valid_campaign_present") is not True
    or receipt.get("test_role_predictions_present") is not False
    or receipt.get("fixed_seed_rectangle_complete_fits") != 4
    or receipt.get("fixed_seed_rectangle_missing_fits") != 11
    or receipt.get("development_gate_passed") is not False
    or receipt.get("global_profile_fixture_admission_complete") is not False
    or receipt.get("observed_multiome_task_eligible") is not False
    or receipt.get("valid_next_centrally_dispatched_gpu_unit_exists") is not False
    or receipt.get("new_gpu_queue_registration_authorized") is not False
    or receipt.get("direct_gpu_submission_authorized") is not False
    or receipt.get("prediction_values_read") is not False
    or receipt.get("metric_values_read") is not False
    or receipt.get("sealed_labels_or_outcomes_read") is not False
    or receipt.get("gpu_queue_registered") is not False
    or receipt.get("gpu_job_submitted") is not False
    or "Ran 5 tests" not in stderr
    or not stderr.rstrip().endswith("OK")
):
    raise SystemExit("scBasset production-readiness receipt differs")

freeze_tree(
    root,
    {
        "artifact_class": "scbasset_production_readiness",
        "model_id": "scbasset",
        "artifact_role": "outcome_blind_semantics_coverage_leakage_and_dispatch_reconciliation",
        "tests": 5,
        "native_observed_atac_role": "outer_training_supervision_only",
        "observed_query_atac_consumed": False,
        "inductive_observed_multiome_model": False,
        "five_outer_fold_valid_campaign_present": True,
        "test_role_predictions_present": False,
        "fixed_seed_rectangle_complete_fits": 4,
        "fixed_seed_rectangle_missing_fits": 11,
        "development_gate_passed": False,
        "global_profile_fixture_admission_complete": False,
        "observed_multiome_task_eligible": False,
        "valid_next_centrally_dispatched_gpu_unit_exists": False,
        "gpu_queue_registered": False,
        "gpu_job_submitted": False,
        "prediction_values_read": False,
        "metric_values_read": False,
        "sealed_labels_or_outcomes_read": False,
        "contract_sha256": receipt["contract_sha256"],
        "status": "passed",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'scBasset production-readiness validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
