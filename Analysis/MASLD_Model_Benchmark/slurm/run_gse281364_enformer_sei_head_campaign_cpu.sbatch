#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-train-604
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
CONFIG=${ROOT}/config/gse281364_enformer_sei_head_campaign.json
TARGET=${ROOT}/executions/model-cpu-train-604-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Enformer/Sei head campaign: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
export PYTHONPATH=${ROOT}:${ROOT}/src
cd "${ROOT}"

sha256sum \
  config/gse281364_enformer_sei_head_campaign.json \
  scripts/freeze_gse281364_enformer_sei_head_taskspec.py \
  scripts/fit_gse281364_enformer_sei_heads.py \
  scripts/evaluate_gse281364_enformer_sei_heads.py \
  tests/unit/test_gse281364_enformer_sei_head_campaign.py \
  slurm/run_gse281364_enformer_sei_head_campaign_cpu.sbatch \
  > "${STAGE}/source.start.sha256"

"${PYTHON}" -m unittest -v \
  tests.unit.test_gse281364_enformer_sei_head_campaign \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

# TaskSpec is frozen before any reporter outcome is opened.
"${PYTHON}" "${ROOT}/scripts/freeze_gse281364_enformer_sei_head_taskspec.py" \
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
    "artifact_class": "gse281364_enformer_sei_head_taskspec",
    "dataset_id": "gse281364",
    "task_id": "variant_to_regulation",
    "candidate_count": 6,
    "outcomes_read": False,
    "model_fit": False,
    "sealed_assets_read": False,
    "status": "pass_prespecified_outcome_blind_taskspec",
})
verify_frozen_tree(root)
PY
TASKSPEC_SHA=$(sha256sum "${STAGE}/task_spec/ARTIFACTS.json" | awk '{print $1}')

"${PYTHON}" "${ROOT}/scripts/fit_gse281364_enformer_sei_heads.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --task-spec-tree "${STAGE}/task_spec" \
  --task-spec-sha256 "${TASKSPEC_SHA}" \
  --output "${STAGE}/fit" \
  > "${STAGE}/validation/fit.stdout.json" \
  2> "${STAGE}/validation/fit.stderr.txt"

"${PYTHON}" - "${STAGE}/fit" "${TASKSPEC_SHA}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(root, {
    "artifact_class": "gse281364_enformer_sei_restricted_head_fit",
    "dataset_id": "gse281364",
    "task_id": "variant_to_regulation",
    "candidate_count": 6,
    "prediction_rows": 61980,
    "task_spec_artifacts_sha256": sys.argv[2],
    "metrics_calculated": False,
    "sealed_assets_read": False,
    "status": "pass_restricted_development_oof_fit",
})
verify_frozen_tree(root)
PY
FIT_SHA=$(sha256sum "${STAGE}/fit/ARTIFACTS.json" | awk '{print $1}')

"${PYTHON}" "${ROOT}/scripts/evaluate_gse281364_enformer_sei_heads.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --task-spec-tree "${STAGE}/task_spec" \
  --task-spec-sha256 "${TASKSPEC_SHA}" \
  --fit-root "${STAGE}/fit" \
  --fit-sha256 "${FIT_SHA}" \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/validation/evaluation.stdout.json" \
  2> "${STAGE}/validation/evaluation.stderr.txt"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  config/gse281364_enformer_sei_head_campaign.json \
  scripts/freeze_gse281364_enformer_sei_head_taskspec.py \
  scripts/fit_gse281364_enformer_sei_heads.py \
  scripts/evaluate_gse281364_enformer_sei_heads.py \
  tests/unit/test_gse281364_enformer_sei_head_campaign.py \
  slurm/run_gse281364_enformer_sei_head_campaign_cpu.sbatch \
  > "${STAGE}/source.terminal.sha256"
cmp --silent "${STAGE}/source.start.sha256" "${STAGE}/source.terminal.sha256"
sha256sum --check --strict "${STAGE}/source.terminal.sha256" \
  > "${STAGE}/validation/terminal_source_recheck.stdout.txt"

"${PYTHON}" - "${STAGE}" "${TASKSPEC_SHA}" "${FIT_SHA}" <<'PY'
import json
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

stage = Path(sys.argv[1]).resolve(strict=True)
fit = json.loads((stage / "fit/receipt.json").read_text(encoding="utf-8"))
evaluation = json.loads((stage / "evaluation/receipt.json").read_text(encoding="utf-8"))
if (
    fit["status"] != "pass_restricted_development_oof_fit"
    or fit["metrics_calculated"]
    or fit["sealed_assets_read"]
    or evaluation["status"] != "pass_independent_restricted_development_evaluation"
    or evaluation["models_ranked_as_interchangeable"]
    or evaluation["shortlist_created"]
    or evaluation["champion_claim"]
    or evaluation["sealed_assets_read"]
):
    raise SystemExit("terminal campaign firewall differs")
freeze_tree(stage, {
    "artifact_class": "gse281364_enformer_sei_restricted_head_campaign",
    "dataset_id": "gse281364",
    "task_id": "variant_to_regulation",
    "models": ["enformer_crested_restricted_port", "sei"],
    "candidate_count": 6,
    "prediction_rows": 61980,
    "elements": 1033,
    "long_range_blocks": 239,
    "fixed_seeds": [1103, 2909, 4721, 6673, 8111],
    "task_spec_artifacts_sha256": sys.argv[2],
    "fit_artifacts_sha256": sys.argv[3],
    "family_native_results_reported_separately": True,
    "mandatory_task_native_baselines_complete": False,
    "models_ranked_as_interchangeable": False,
    "shortlist_created": False,
    "champion_claim": False,
    "native_equivalence_claim": False,
    "external_evaluation": False,
    "sealed_assets_read": False,
    "global_frozen_census_changed": False,
    "status": "passed_restricted_development_head_campaign",
})
verify_frozen_tree(stage)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Completed Enformer/Sei development head campaign: %s\n' "${TARGET}"
