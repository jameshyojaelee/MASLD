#!/usr/bin/env bash
#SBATCH --job-name=model-work-502
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
: "${VALIDATION_ROOT:?submit with the exact frozen v2 validation root}"
: "${VALIDATION_SHA:?submit with the exact frozen v2 validation ARTIFACTS SHA-256}"
export VALIDATION_ROOT
export VALIDATION_SHA

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
CONFIG=${ROOT}/config/gse281364_open_sequence_head_campaign.toml
TARGET=${ROOT}/executions/model-work-502-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE281364 head campaign v2: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
export PYTHONPATH=${ROOT}:${ROOT}/src
cd "${ROOT}"

# This exact frozen-tree and source check precedes every TaskSpec authority check,
# unit test, outcome read, model fit, and prediction write in production v2.
"${PYTHON}" "${ROOT}/scripts/verify_gse281364_open_sequence_production_binding.py" \
  --root "${ROOT}" \
  --validation-root "${VALIDATION_ROOT}" \
  --validation-sha "${VALIDATION_SHA}" \
  > "${STAGE}/validation/binding.stdout.json" \
  2> "${STAGE}/validation/binding.stderr.txt"
cp "${VALIDATION_ROOT}/source.sha256" "${STAGE}/validation/validated_source.sha256"
printf '%s  %s\n' "${VALIDATION_SHA}" "${VALIDATION_ROOT}/ARTIFACTS.json" \
  > "${STAGE}/validation/validation_artifacts.sha256"

"${PYTHON}" -m unittest -v \
  tests.unit.test_gse281364_open_sequence_head_campaign \
  tests.unit.test_verify_gse281364_open_sequence_production_binding \
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
  config/gse281364_open_sequence_head_campaign.toml \
  scripts/freeze_gse281364_open_sequence_taskspec.py \
  scripts/fit_gse281364_open_sequence_heads.py \
  scripts/evaluate_gse281364_open_sequence_heads.py \
  scripts/verify_gse281364_open_sequence_production_binding.py \
  tests/unit/test_gse281364_open_sequence_head_campaign.py \
  tests/unit/test_verify_gse281364_open_sequence_production_binding.py \
  slurm/run_gse281364_open_sequence_head_campaign_v2_cpu.sbatch \
  slurm/validate_gse281364_open_sequence_head_campaign_v2_cpu.sbatch \
  > "${STAGE}/source.sha256"

# Detect any config, code, test, or wrapper drift that occurred after startup.
cmp --silent "${STAGE}/source.sha256" "${VALIDATION_ROOT}/source.sha256"
sha256sum --check --strict "${STAGE}/source.sha256" \
  > "${STAGE}/validation/terminal_source_recheck.stdout.txt"

"${PYTHON}" - "${STAGE}" "${TASKSPEC_SHA}" "${VALIDATION_ROOT}" "${VALIDATION_SHA}" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


stage = Path(sys.argv[1]).resolve(strict=True)
task_spec_sha256 = sys.argv[2]
validation_root = Path(sys.argv[3]).resolve(strict=True)
validation_sha256 = sys.argv[4]
binding = json.loads(
    (stage / "validation/binding.stdout.json").read_text(encoding="utf-8")
)
fit = json.loads((stage / "fit/receipt.json").read_text(encoding="utf-8"))
evaluation = json.loads((stage / "evaluation/receipt.json").read_text(encoding="utf-8"))
terminal_checks = (
    stage / "validation/terminal_source_recheck.stdout.txt"
).read_text(encoding="utf-8").splitlines()
if (
    binding["status"] != "pass_exact_frozen_validation_binding"
    or binding["validation_root"] != validation_root.as_posix()
    or binding["validation_artifacts_sha256"] != validation_sha256
    or binding["outcomes_read"]
    or binding["mandatory_baselines_complete"]
    or not binding["shortlist_blocked"]
    or not binding["complementarity_blocked"]
    or len(terminal_checks) != 9
    or not all(line.endswith(": OK") for line in terminal_checks)
    or fit["status"] != "pass_prespecified_exposed_development_oof_fit"
    or fit["prediction_rows"] != 123960
    or fit["task_spec_artifacts_sha256"] != task_spec_sha256
    or fit["biological_donors"] != 0
    or fit["experimental_replicates_used_as_independent_donors"]
    or fit["metrics_calculated"]
    or fit["sealed_assets_read"]
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
    or evaluation["mandatory_baselines_complete"]
    or len(evaluation["mandatory_baselines_missing"]) != 5
    or not evaluation["shortlist_blocked"]
    or evaluation["development_shortlist_count"] != 0
    or not evaluation["complementarity_blocked"]
    or not evaluation["conditional_trigger_blocked"]
):
    raise SystemExit("GSE281364 head-campaign v2 production receipt differs")
freeze_tree(stage, {
    "artifact_class": "gse281364_open_sequence_head_campaign_v2",
    "dataset_id": "gse281364",
    "task_id": "variant_to_regulation",
    "task_spec_artifacts_sha256": task_spec_sha256,
    "validation_root": validation_root.as_posix(),
    "validation_artifacts_sha256": validation_sha256,
    "validation_source_manifest_sha256": binding["source_manifest_sha256"],
    "production_wrapper_sha256": binding["production_wrapper_sha256"],
    "terminal_source_recheck_passed": True,
    "candidate_count": 12,
    "prediction_rows": 123960,
    "biological_donors": 0,
    "experimental_replicates_per_context": 4,
    "claim_role": "exposed_development_MPRA_proxy_only",
    "campaign_scope": "partial_common_head_campaign",
    "mandatory_baselines_complete": False,
    "mandatory_baselines_missing": binding["mandatory_baselines_missing"],
    "shortlist_blocked": True,
    "finalist_claim_blocked": True,
    "complementarity_blocked": True,
    "conditional_trigger_blocked": True,
    "champion_claim": False,
    "external_evaluation": False,
    "sealed_assets_read": False,
    "status": "pass_independent_exposed_development_evaluation",
})
verify_frozen_tree(stage)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 open-sequence head campaign v2: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
