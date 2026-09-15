#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-score-511
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${VIEWS_ROOT:?VIEWS_ROOT is required}"
: "${VIEWS_SHA256:?VIEWS_SHA256 is required}"
: "${PREDICTIONS_ROOT:?PREDICTIONS_ROOT is required}"
: "${PREDICTIONS_SHA256:?PREDICTIONS_SHA256 is required}"

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
TASK=${ROOT}/config/evaluation/paired_bulk_rna_h3k27ac_task.toml
TASK_SHA256=8bf05fe693c787ffdd8fe11443e019e668ecd03318c708e390785edd837b1a28
GATE=${ROOT}/config/evaluation/paired_bulk_rna_h3k27ac_promotion_gate.json
GATE_SHA256=f26a2048346cb60f810407971c41bbc22fe862444ee5438c90da954eb44b21a1
TARGET=${ROOT}/executions/model-cpu-score-511-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite paired baseline evaluation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${VIEWS_SHA256}" "${VIEWS_ROOT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${PREDICTIONS_SHA256}" "${PREDICTIONS_ROOT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${TASK_SHA256}" "${TASK}" | sha256sum --check --strict
printf '%s  %s\n' "${GATE_SHA256}" "${GATE}" | sha256sum --check --strict
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict

cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
"${PYTHON}" -m unittest -v tests.unit.test_gse267145_paired_baseline_campaign \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${PYTHON}" - "${VIEWS_ROOT}" "${PREDICTIONS_ROOT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
PY
EVALUATOR_VIEW=${VIEWS_ROOT}/evaluator
EVALUATOR_SHA256=$(sha256sum "${EVALUATOR_VIEW}/ARTIFACTS.json" | awk '{print $1}')
"${PYTHON}" scripts/evaluate_gse267145_paired_baseline_predictions.py \
  --evaluator-view "${EVALUATOR_VIEW}" \
  --evaluator-artifacts-sha256 "${EVALUATOR_SHA256}" \
  --predictions "${PREDICTIONS_ROOT}" \
  --prediction-artifacts-sha256 "${PREDICTIONS_SHA256}" \
  --task-spec "${TASK}" \
  --task-spec-sha256 "${TASK_SHA256}" \
  --promotion-gate "${GATE}" \
  --promotion-gate-sha256 "${GATE_SHA256}" \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/validation/evaluator_stdout.json" \
  2> "${STAGE}/validation/evaluator_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/evaluate_gse267145_paired_baseline_predictions.py \
  tests/unit/test_gse267145_paired_baseline_campaign.py \
  slurm/evaluate_gse267145_paired_baselines_cpu.sbatch \
  "${VIEWS_ROOT}/ARTIFACTS.json" \
  "${EVALUATOR_VIEW}/ARTIFACTS.json" \
  "${PREDICTIONS_ROOT}/ARTIFACTS.json" \
  "${TASK}" \
  "${GATE}" \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
verify_frozen_tree(root / "evaluation")
receipt = json.loads((root / "evaluation/evaluation_receipt.json").read_text())
if (
    receipt["status"] != "pass_development_evaluation_only"
    or receipt["participants"] != 99
    or receipt["biological_unit"] != "participant"
    or receipt["models_fit_inside_evaluator"]
    or receipt["preprocessing_fit_inside_evaluator"]
    or not receipt["held_participant_h3_read_only_inside_evaluator"]
    or receipt["histology_sex_age_nas_fibrosis_read"]
    or receipt["coordinate_interpretation_used"]
    or receipt["model_selection_performed"]
    or receipt["champion_promotion_performed"]
    or receipt["external_evaluation_performed"]
    or receipt["bootstrap"]["replicates"] != 10_000
):
    raise SystemExit("paired baseline evaluation receipt differs")
freeze_tree(root, {
    "artifact_class": "gse267145_paired_baseline_evaluation_campaign",
    "participants": 99,
    "models": 5,
    "model_selection_performed": False,
    "champion_promotion_performed": False,
    "external_evaluation_performed": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/evaluation/ARTIFACTS.json"
