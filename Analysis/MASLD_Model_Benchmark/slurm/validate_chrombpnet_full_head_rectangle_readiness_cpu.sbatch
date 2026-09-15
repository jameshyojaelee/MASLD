#!/usr/bin/env bash
#SBATCH --job-name=model-check-991
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
CONTRACT="${ROOT}/config/artifacts/models/chrombpnet/full_head_rectangle_execution_readiness_20260825.json"
TARGET="${ROOT}/executions/model-check-991-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite ChromBPNet lane validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"

for pair in \
  "b9557709843cfc11b9584957187c7c0d061d950b976315a8ce095e3a0a7dffd3 config/artifacts/models/chrombpnet/full_head_rectangle_execution_readiness_20260825.json" \
  "a44bb212f7bb8a6a4bfc31c18b805a2c7327d9e79d5d91bd7af967c5fa1becf7 scripts/audit_chrombpnet_full_head_rectangle_readiness.py" \
  "53cd2601eb799ee919ec06dd94fdf70a5c56ce4130e60277ad41f59bd836df32 tests/unit/test_audit_chrombpnet_full_head_rectangle_readiness.py"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.unit.test_audit_chrombpnet_full_head_rectangle_readiness \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}/src:${ROOT}" python \
  scripts/audit_chrombpnet_full_head_rectangle_readiness.py \
  --project-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --output "${STAGE}/readiness/receipt.json" \
  > "${STAGE}/validation/audit.stdout.txt"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  config/artifacts/models/chrombpnet/full_head_rectangle_execution_readiness_20260825.json \
  scripts/audit_chrombpnet_full_head_rectangle_readiness.py \
  tests/unit/test_audit_chrombpnet_full_head_rectangle_readiness.py \
  slurm/validate_chrombpnet_full_head_rectangle_readiness_cpu.sbatch \
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
    != "pass_chrombpnet_full_head_lane_existing_central_queue_authorized"
    or receipt.get("one_fit_two_views") is not True
    or receipt.get("primary_biological_view") != "chrombpnet_nobias"
    or receipt.get("assay_qc_view") != "chrombpnet_full"
    or receipt.get("expected_fits") != 125
    or receipt.get("complete_fits_at_freeze") != 19
    or receipt.get("missing_fits_at_freeze") != 106
    or receipt.get("complete_prediction_views_at_freeze") != 38
    or receipt.get("missing_prediction_views_at_freeze") != 212
    or receipt.get("full_compatible_rectangle") is not False
    or receipt.get("missing_fits_already_centrally_registered") != 106
    or receipt.get("central_task_bundles") != 7
    or receipt.get("new_gpu_queue_registration_needed") is not False
    or receipt.get("direct_gpu_submission_authorized") is not False
    or receipt.get("gpu_squeue_ceiling") != 5
    or receipt.get("maximum_running_gpu_jobs") != 4
    or receipt.get("maximum_pending_gpu_jobs") != 1
    or receipt.get("qos") != "nslab"
    or receipt.get("license_authority_contradiction_present") is not True
    or receipt.get("weight_redistribution_authorized") is not False
    or receipt.get("open_champion_eligible_now") is not False
    or receipt.get("partial_results_ranked") is not False
    or receipt.get("prediction_values_read") is not False
    or receipt.get("metric_values_read") is not False
    or receipt.get("development_outcomes_read") is not False
    or receipt.get("sealed_labels_or_outcomes_read") is not False
    or receipt.get("gpu_job_submitted") is not False
    or "Ran 5 tests" not in stderr
    or not stderr.rstrip().endswith("OK")
):
    raise SystemExit("ChromBPNet full-head lane readiness receipt differs")

freeze_tree(
    root,
    {
        "artifact_class": "chrombpnet_full_head_rectangle_execution_readiness",
        "model_id": "chrombpnet",
        "lane_id": "chrombpnet_full_head_two_view_five_seed_rectangle",
        "artifact_role": "outcome_blind_two_view_rectangle_and_dispatch_coverage",
        "tests": 5,
        "one_fit_two_views": True,
        "primary_biological_view": "chrombpnet_nobias",
        "assay_qc_view": "chrombpnet_full",
        "expected_fits": 125,
        "complete_fits_at_freeze": 19,
        "missing_fits_at_freeze": 106,
        "complete_prediction_views_at_freeze": 38,
        "missing_prediction_views_at_freeze": 212,
        "full_compatible_rectangle": False,
        "missing_fits_already_centrally_registered": 106,
        "central_task_bundles": 7,
        "new_gpu_queue_registration_needed": False,
        "gpu_job_submitted": False,
        "gpu_squeue_ceiling": 5,
        "maximum_running_gpu_jobs": 4,
        "maximum_pending_gpu_jobs": 1,
        "qos": "nslab",
        "license_authority_contradiction_present": True,
        "weight_redistribution_authorized": False,
        "open_champion_eligible_now": False,
        "partial_results_ranked": False,
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
printf 'ChromBPNet full-head rectangle readiness: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
