#!/usr/bin/env bash
#SBATCH --job-name=model-check-218
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
CONFIG=${ROOT}/config/local_variant_family_reconciliation.json
TARGET=${ROOT}/executions/model-check-218-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite local variant reconciliation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
export PYTHONPATH=${ROOT}:${ROOT}/src

"${PYTHON}" -m unittest -v \
  tests.unit.test_preflight_local_variant_family_reconciliation \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

"${PYTHON}" "${ROOT}/scripts/preflight_local_variant_family_reconciliation.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/reconciliation" \
  > "${STAGE}/validation/preflight.stdout.json" \
  2> "${STAGE}/validation/preflight.stderr.txt"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/preflight_local_variant_family_reconciliation.py" \
  "${ROOT}/tests/unit/test_preflight_local_variant_family_reconciliation.py" \
  "${ROOT}/slurm/preflight_local_variant_family_reconciliation_cpu.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree


stage = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads(
    (stage / "reconciliation/receipt.json").read_text(encoding="utf-8")
)
if (
    receipt["status"] != "pass_reconciled_native_tasks_no_model_promoted"
    or receipt["model_rows"] != 10
    or receipt["native_task_rows"] != 5
    or receipt["authority_rows"] != 38
    or receipt["additional_development_screens_allowed_now"]
    or len(receipt["conditional_next_screens"]) != 3
    or receipt["active_campaign_outputs_bound_or_read"]
    or receipt["sealed_assets_read"]
    or receipt["test_outcomes_read"]
    or receipt["new_training_scoring_or_model_import_performed"]
    or receipt["primary_signed_eQTL_or_ieQTL_candidate_present"]
    or receipt["champion_claim_allowed"]
    or receipt["universal_model_claim_allowed"]
    or receipt["clinical_claim_allowed"]
):
    raise SystemExit("local variant family reconciliation receipt differs")

freeze_tree(
    stage,
    {
        "artifact_class": "local_variant_family_reconciliation",
        "model_count": 10,
        "native_task_count": 5,
        "active_campaign_outputs_bound_or_read": False,
        "sealed_assets_read": False,
        "test_outcomes_read": False,
        "new_training_scoring_or_model_import_performed": False,
        "champion_claim_allowed": False,
        "status": "pass_reconciled_native_tasks_no_model_promoted",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Local variant family reconciliation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
