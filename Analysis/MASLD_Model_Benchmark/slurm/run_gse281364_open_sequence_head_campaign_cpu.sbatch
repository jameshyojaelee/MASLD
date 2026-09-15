#!/usr/bin/env bash
#SBATCH --job-name=model-work-501
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=36:00:00
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
TARGET=${ROOT}/executions/model-work-501-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE281364 head campaign: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
export PYTHONPATH=${ROOT}:${ROOT}/src

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
TASKSPEC_SHA=$(sha256sum "${STAGE}/task_spec/ARTIFACTS.json" | awk '{print $1}')

"${PYTHON}" "${ROOT}/scripts/fit_gse281364_open_sequence_heads.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --task-spec-tree "${STAGE}/task_spec" \
  --task-spec-sha256 "${TASKSPEC_SHA}" \
  --output "${STAGE}/fit" \
  > "${STAGE}/validation/fit.stdout.json" \
  2> "${STAGE}/validation/fit.stderr.txt"

"${PYTHON}" "${ROOT}/scripts/evaluate_gse281364_open_sequence_heads.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --task-spec-sha256 "${TASKSPEC_SHA}" \
  --fit-root "${STAGE}/fit" \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/validation/evaluation.stdout.json" \
  2> "${STAGE}/validation/evaluation.stderr.txt"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/freeze_gse281364_open_sequence_taskspec.py" \
  "${ROOT}/scripts/fit_gse281364_open_sequence_heads.py" \
  "${ROOT}/scripts/evaluate_gse281364_open_sequence_heads.py" \
  "${ROOT}/tests/unit/test_gse281364_open_sequence_head_campaign.py" \
  "${ROOT}/slurm/run_gse281364_open_sequence_head_campaign_cpu.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" "${TASKSPEC_SHA}" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


stage = Path(sys.argv[1]).resolve(strict=True)
task_spec_sha256 = sys.argv[2]
fit = json.loads((stage / "fit/receipt.json").read_text(encoding="utf-8"))
evaluation = json.loads((stage / "evaluation/receipt.json").read_text(encoding="utf-8"))
if (
    fit["status"] != "pass_prespecified_exposed_development_oof_fit"
    or fit["prediction_rows"] != 123960
    or fit["task_spec_artifacts_sha256"] != task_spec_sha256
    or fit["biological_donors"] != 0
    or fit["experimental_replicates_used_as_independent_donors"]
    or fit["metrics_calculated"]
    or fit["sealed_assets_read"]
    or fit["stack_fit"]
    or fit["residual_correlation_calculated"]
    or fit["conditional_model_built_or_fit"]
    or fit["mandatory_baselines_complete"]
    or len(fit["mandatory_baselines_missing"]) != 5
    or not fit["shortlist_blocked"]
    or not fit["complementarity_blocked"]
    or evaluation["status"] != "pass_independent_exposed_development_evaluation"
    or evaluation["prediction_rows"] != 123960
    or evaluation["biological_donors"] != 0
    or evaluation["experimental_replicates_used_as_independent_donors"]
    or evaluation["confirmatory_p_values_reported"]
    or evaluation["champion_claim"]
    or evaluation["external_evaluation"]
    or evaluation["sealed_assets_read"]
    or evaluation["stack_fit"]
    or evaluation["residual_correlation_calculated"]
    or evaluation["conditional_model_built_or_fit"]
    or evaluation["mandatory_baselines_complete"]
    or len(evaluation["mandatory_baselines_missing"]) != 5
    or not evaluation["shortlist_blocked"]
    or evaluation["development_shortlist_count"] != 0
    or not evaluation["complementarity_blocked"]
    or not evaluation["conditional_trigger_blocked"]
):
    raise SystemExit("GSE281364 head-campaign production receipt differs")
freeze_tree(stage, {
    "artifact_class": "gse281364_open_sequence_head_campaign",
    "dataset_id": "gse281364",
    "task_id": "variant_to_regulation",
    "task_spec_artifacts_sha256": task_spec_sha256,
    "candidate_count": 12,
    "prediction_rows": 123960,
    "biological_donors": 0,
    "experimental_replicates_per_context": 4,
    "claim_role": "exposed_development_MPRA_proxy_only",
    "campaign_scope": "partial_common_head_campaign",
    "mandatory_baselines_complete": False,
    "shortlist_blocked": True,
    "complementarity_blocked": True,
    "champion_claim": False,
    "external_evaluation": False,
    "sealed_assets_read": False,
    "status": "pass_independent_exposed_development_evaluation",
})
verify_frozen_tree(stage)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 open-sequence head campaign: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
