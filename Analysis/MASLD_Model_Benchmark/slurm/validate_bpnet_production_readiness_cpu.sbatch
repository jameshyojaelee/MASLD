#!/usr/bin/env bash
#SBATCH --job-name=model-check-990
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
CONTRACT="${ROOT}/config/artifacts/models/bpnet/production_readiness_20260825.json"
TARGET="${ROOT}/executions/model-check-990-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite BPNet readiness validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"

for pair in \
  "00fd03eb0fd7ca1dead56edad457c60fa447a7c5a84def8e11dc9bf7218cddd8 config/artifacts/models/bpnet/production_readiness_20260825.json" \
  "034700f5a060acb7f4b2ca40a6f05fe57c1dedda968452b43c2c16a064385811 scripts/audit_bpnet_production_readiness.py" \
  "d052ca895fa7ef0b7cecf93697d9c219105de99984c109f133b160cbe9a6c349 tests/unit/test_audit_bpnet_production_readiness.py"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.unit.test_audit_bpnet_production_readiness \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}/src:${ROOT}" python \
  scripts/audit_bpnet_production_readiness.py \
  --project-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --output "${STAGE}/readiness/receipt.json" \
  > "${STAGE}/validation/audit.stdout.txt"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  config/artifacts/models/bpnet/production_readiness_20260825.json \
  scripts/audit_bpnet_production_readiness.py \
  tests/unit/test_audit_bpnet_production_readiness.py \
  slurm/validate_bpnet_production_readiness_cpu.sbatch \
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
    receipt.get("status") != "pass_bpnet_reconciled_existing_central_queue_authorized"
    or receipt.get("bpnet_is_separate_architecture_from_chrombpnet") is not True
    or receipt.get("bpnet_is_chrombpnet_nobias_view") is not False
    or receipt.get("expected_fits") != 125
    or receipt.get("complete_fits_at_freeze") != 20
    or receipt.get("missing_fits_at_freeze") != 105
    or receipt.get("full_compatible_rectangle") is not False
    or receipt.get("missing_fits_already_centrally_registered") != 105
    or receipt.get("central_task_bundles") != 7
    or receipt.get("new_gpu_queue_registration_needed") is not False
    or receipt.get("direct_gpu_submission_authorized") is not False
    or receipt.get("gpu_squeue_ceiling") != 5
    or receipt.get("maximum_running_gpu_jobs") != 4
    or receipt.get("maximum_pending_gpu_jobs") != 1
    or receipt.get("qos") != "nslab"
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
    raise SystemExit("BPNet production-readiness receipt differs")

freeze_tree(
    root,
    {
        "artifact_class": "bpnet_production_readiness",
        "model_id": "bpnet",
        "artifact_role": "outcome_blind_architecture_rectangle_license_and_dispatch_reconciliation",
        "tests": 5,
        "separate_from_chrombpnet_and_nobias_view": True,
        "expected_fits": 125,
        "complete_fits_at_freeze": 20,
        "missing_fits_at_freeze": 105,
        "full_compatible_rectangle": False,
        "missing_fits_already_centrally_registered": 105,
        "central_task_bundles": 7,
        "new_gpu_queue_registration_needed": False,
        "gpu_job_submitted": False,
        "gpu_squeue_ceiling": 5,
        "maximum_running_gpu_jobs": 4,
        "maximum_pending_gpu_jobs": 1,
        "qos": "nslab",
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
printf 'BPNet production-readiness validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
