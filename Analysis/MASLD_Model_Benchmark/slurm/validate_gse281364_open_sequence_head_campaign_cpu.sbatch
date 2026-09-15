#!/usr/bin/env bash
#SBATCH --job-name=model-check-223
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
CONFIG=${ROOT}/config/gse281364_open_sequence_head_campaign.toml
TARGET=${ROOT}/executions/model-check-223-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE281364 head-campaign validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
export PYTHONPATH=${ROOT}:${ROOT}/src

"${PYTHON}" -m py_compile \
  "${ROOT}/scripts/freeze_gse281364_open_sequence_taskspec.py" \
  "${ROOT}/scripts/fit_gse281364_open_sequence_heads.py" \
  "${ROOT}/scripts/evaluate_gse281364_open_sequence_heads.py" \
  "${ROOT}/tests/unit/test_gse281364_open_sequence_head_campaign.py"

"${PYTHON}" -m unittest -v \
  tests.unit.test_gse281364_open_sequence_head_campaign \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

"${PYTHON}" "${ROOT}/scripts/freeze_gse281364_open_sequence_taskspec.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/task_spec" \
  > "${STAGE}/validation/taskspec.stdout.json" \
  2> "${STAGE}/validation/taskspec.stderr.txt"

"${PYTHON}" - "${STAGE}/task_spec" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(root, {
    "artifact_class": "gse281364_open_sequence_taskspec",
    "dataset_id": "gse281364",
    "task_id": "variant_to_regulation",
    "candidate_count": 12,
    "biological_donors": 0,
    "experimental_replicates_per_context": 4,
    "outcomes_read": False,
    "model_fit": False,
    "status": "pass_prespecified_taskspec",
})
verify_frozen_tree(root)
PY

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/freeze_gse281364_open_sequence_taskspec.py" \
  "${ROOT}/scripts/fit_gse281364_open_sequence_heads.py" \
  "${ROOT}/scripts/evaluate_gse281364_open_sequence_heads.py" \
  "${ROOT}/tests/unit/test_gse281364_open_sequence_head_campaign.py" \
  "${ROOT}/slurm/validate_gse281364_open_sequence_head_campaign_cpu.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


stage = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((stage / "task_spec/receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_prespecified_taskspec"
    or receipt["candidate_count"] != 12
    or receipt["biological_donors"] != 0
    or receipt["measurement_replicates_per_context"] != 4
    or receipt["measurement_replicates_are_biological_replicates"]
    or receipt["outcome_values_read"]
    or receipt["model_fit"]
    or receipt["predictions_generated"]
    or receipt["stack_fit"]
    or receipt["residual_correlation_calculated"]
    or receipt["conditional_model_built_or_fit"]
    or receipt["sealed_assets_read"]
    or receipt["mandatory_baselines_complete"]
    or len(receipt["mandatory_baselines_missing"]) != 5
    or not receipt["shortlist_blocked"]
    or not receipt["finalist_claim_blocked"]
    or not receipt["complementarity_blocked"]
    or not receipt["conditional_trigger_blocked"]
):
    raise SystemExit("GSE281364 head-campaign validation receipt differs")
freeze_tree(stage, {
    "artifact_class": "gse281364_open_sequence_head_campaign_validation",
    "dataset_id": "gse281364",
    "task_id": "variant_to_regulation",
    "unit_tests": 7,
    "candidate_count": 12,
    "biological_donors": 0,
    "outcomes_read": False,
    "model_fit": False,
    "predictions_generated": False,
    "campaign_scope": "partial_common_head_campaign",
    "mandatory_baselines_complete": False,
    "shortlist_blocked": True,
    "complementarity_blocked": True,
    "sealed_assets_read": False,
    "status": "passed",
})
verify_frozen_tree(stage)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 open-sequence head validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
