#!/usr/bin/env bash
#SBATCH --job-name=model-check-219
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=04:00:00
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
CONFIG=${ROOT}/config/complementarity_readiness_audit.toml
TARGET=${ROOT}/executions/model-check-219-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite complementarity readiness audit: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
export PYTHONPATH=${ROOT}:${ROOT}/src

"${PYTHON}" -m unittest -v \
  tests.unit.test_audit_complementarity_readiness \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

"${PYTHON}" "${ROOT}/scripts/audit_complementarity_readiness.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/audit" \
  > "${STAGE}/validation/audit.stdout.json" \
  2> "${STAGE}/validation/audit.stderr.txt"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/audit_complementarity_readiness.py" \
  "${ROOT}/tests/unit/test_audit_complementarity_readiness.py" \
  "${ROOT}/slurm/audit_complementarity_readiness_cpu.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree


stage = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((stage / "audit/receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "blocked_no_trigger_compatible_prediction_pair"
    or receipt["authority_rows"] != 19
    or receipt["prediction_lane_rows"] != 8
    or receipt["pair_audit_rows"] != 6
    or receipt["required_stacking_source_rows"] != 6
    or receipt["bound_stacking_source_rows"] != 0
    or receipt["trigger_ready_pair_rows"] != 0
    or len(receipt["scbasset_midas_row_overlap"]) != 5
    or receipt["residual_correlation_calculated"]
    or receipt["stack_weights_fit"]
    or receipt["prediction_values_read"]
    or not receipt["row_identifiers_read"]
    or receipt["raw_or_summary_outcomes_read"]
    or receipt["sealed_assets_read"]
    or receipt["novel_model_built_or_fit"]
    or receipt["conditional_model_build_authorized"]
):
    raise SystemExit("complementarity readiness receipt differs")

freeze_tree(
    stage,
    {
        "artifact_class": "complementarity_readiness_audit",
        "prediction_lane_count": 8,
        "pair_audit_count": 6,
        "bound_stacking_source_count": 0,
        "trigger_ready_pair_count": 0,
        "prediction_values_read": False,
        "row_identifiers_read": True,
        "raw_or_summary_outcomes_read": False,
        "sealed_assets_read": False,
        "novel_model_built_or_fit": False,
        "conditional_model_build_authorized": False,
        "status": "blocked_no_trigger_compatible_prediction_pair",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Complementarity readiness audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
