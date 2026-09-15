#!/usr/bin/env bash
#SBATCH --job-name=model-check-987
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
CONTRACT="${ROOT}/config/artifacts/models/regulatory_native_controls/production_readiness_20260825.json"
TARGET="${ROOT}/executions/model-check-987-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite sequence-control readiness validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"

for pair in \
  "78603ed859f1247c213108fcfa1c281d464c4f0871dc135ee7d24db60ab58f37 config/artifacts/models/regulatory_native_controls/production_readiness_20260825.json" \
  "3de616a6150897a7c49a214abfdc9a866b04fd2d0a6c57d347525f98a994c08d scripts/audit_sequence_control_production_readiness.py" \
  "09318a5bbcb4b075c8b46ef33d9f154b13a588032aa38e3623f467c3a3a74415 tests/unit/test_audit_sequence_control_production_readiness.py"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.unit.test_audit_sequence_control_production_readiness \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}/src:${ROOT}" python \
  scripts/audit_sequence_control_production_readiness.py \
  --project-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --output "${STAGE}/readiness/receipt.json" \
  > "${STAGE}/validation/audit.stdout.txt"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  config/artifacts/models/regulatory_native_controls/production_readiness_20260825.json \
  scripts/audit_sequence_control_production_readiness.py \
  tests/unit/test_audit_sequence_control_production_readiness.py \
  slurm/validate_sequence_control_production_readiness_cpu.sbatch \
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
    receipt.get("status")
    != "pass_controls_reconciled_waiting_for_full_rectangles_license_and_external_evaluation"
    or receipt.get("controls_are_separate_architectures") is not True
    or receipt.get("controls_share_checkpoint") is not False
    or receipt.get("parameter_budget_matched") is not True
    or receipt.get("cnn_complete_fits_at_freeze") != 12
    or receipt.get("cnn_missing_fits_at_freeze") != 113
    or receipt.get("transformer_complete_fits_at_freeze") != 11
    or receipt.get("transformer_missing_fits_at_freeze") != 114
    or receipt.get("full_compatible_rectangle") is not False
    or receipt.get("partial_evaluator_access_authorized") is not False
    or receipt.get("project_repository_license_resolved") is not False
    or receipt.get("open_champion_eligible_now") is not False
    or receipt.get("cnn_missing_fits_already_centrally_registered") != 113
    or receipt.get("transformer_missing_fits_already_centrally_registered") != 114
    or receipt.get("new_gpu_queue_registration_needed") is not False
    or receipt.get("direct_gpu_submission_authorized") is not False
    or receipt.get("gse281364_fits_fill_gse296875_rectangle") is not False
    or receipt.get("prediction_values_read") is not False
    or receipt.get("metric_values_read") is not False
    or receipt.get("development_outcomes_read") is not False
    or receipt.get("sealed_labels_or_outcomes_read") is not False
    or receipt.get("partial_results_ranked") is not False
    or receipt.get("gpu_job_submitted") is not False
    or "Ran 6 tests" not in stderr
    or not stderr.rstrip().endswith("OK")
):
    raise SystemExit("sequence-control production-readiness receipt differs")

freeze_tree(
    root,
    {
        "artifact_class": "sequence_control_production_readiness",
        "model_ids": ["sequence_cnn_control", "sequence_transformer_control"],
        "artifact_role": "outcome_blind_architecture_rectangle_split_license_and_dispatch_reconciliation",
        "tests": 6,
        "controls_are_separate_architectures": True,
        "controls_share_checkpoint": False,
        "parameter_budget_matched": True,
        "cnn_complete_fits_at_freeze": 12,
        "cnn_missing_fits_at_freeze": 113,
        "transformer_complete_fits_at_freeze": 11,
        "transformer_missing_fits_at_freeze": 114,
        "full_compatible_rectangle": False,
        "partial_evaluator_access_authorized": False,
        "project_repository_license_resolved": False,
        "open_champion_eligible_now": False,
        "missing_control_fits_already_centrally_registered": 227,
        "new_gpu_queue_registration_needed": False,
        "gpu_job_submitted": False,
        "prediction_values_read": False,
        "metric_values_read": False,
        "development_outcomes_read": False,
        "sealed_labels_or_outcomes_read": False,
        "contract_sha256": receipt["contract_sha256"],
        "status": "passed",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Sequence-control production-readiness validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
