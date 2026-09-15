#!/usr/bin/env bash
#SBATCH --job-name=model-eval
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-eval-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-eval-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

: "${COMPLETION_ROOT:?completion root must be exported by the exact-coverage gate}"
: "${COMPLETION_SHA:?completion ARTIFACTS hash must be exported by the exact-coverage gate}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PREPARE=${ROOT}/scripts/prepare_sequence_evaluation_matrix.py
BASE_EVALUATOR=${ROOT}/scripts/evaluate_sequence_accessibility.py
EVALUATOR=${ROOT}/scripts/evaluate_sequence_task_native_five_seed_rectangle.py
ACCEPTANCE=${ROOT}/executions/sequence-regulatory-completion-gate-tests-21100975
TARGET=${ROOT}/executions/sequence-task-native-five-seed-evaluation-${SLURM_JOB_ID}
STAGE=${TARGET}.staging

preserve_failed_stage() {
  local failed=${TARGET}.failed
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite task-native evaluation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

cd "${ROOT}"
printf '%s  %s\n' \
  5b835e722619e232db0aa0731cfededaa4494a4e117b0e609cf816f082a82258 "${PREPARE}" \
  be1b0dac61123d4875aea210d304ae6c2b67605d74079d490fc21f4a07419824 "${BASE_EVALUATOR}" \
  4ce7a88a08c1da6c316b85d8a15b32e637cea4243d0ffadac8883fef79195b12 "${EVALUATOR}" \
  e0b9b25857263a05d855c11a3cbd90dd457acb7876b9f695dde274a5a1252be4 tests/unit/test_prepare_sequence_evaluation_matrix.py \
  0f04f22378eed5a794e7546889bd5d1d2443fb36a3294aac8a6db69e4d3ec589 tests/unit/test_evaluate_sequence_task_native_five_seed_rectangle.py \
  4632e8d6debea40b832365c31fffd301bcd9896dd765f2722f3b2332ff1e4d29 "${ACCEPTANCE}/ARTIFACTS.json" \
  "${COMPLETION_SHA}" "${COMPLETION_ROOT}/ARTIFACTS.json" \
  | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}/src:${ROOT} python - "${ACCEPTANCE}" "${COMPLETION_ROOT}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import verify_frozen_tree

acceptance, completion = map(lambda value: Path(value).resolve(strict=True), sys.argv[1:])
verify_frozen_tree(acceptance)
verify_frozen_tree(completion)
receipt = json.loads(
    (completion / "completion/completion_lock.json").read_text(encoding="utf-8")
)
if (
    receipt["status"] != "ready_both_rectangles_locked"
    or not receipt["evaluation_launch_authorized"]
    or receipt["partial_ranking_authorized"]
    or receipt["prediction_values_read"]
    or receipt["observed_outcomes_read"]
):
    raise SystemExit("exact-coverage completion lock differs")
PY
PYTHONPATH=${ROOT}/src:${ROOT} python "${PREPARE}" \
  --source-manifest "${COMPLETION_ROOT}/completion/sequence_evaluation_sources.tsv" \
  --project-root "${ROOT}" \
  --output "${STAGE}/matrix" \
  > "${STAGE}/validation/prepare.stdout.json" \
  2> "${STAGE}/validation/prepare.stderr.txt"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
PYTHONPATH=${ROOT}/src:${ROOT} python "${EVALUATOR}" \
  --matrix-root "${STAGE}/matrix" \
  --completion-lock "${COMPLETION_ROOT}" \
  --completion-artifacts-sha256 "${COMPLETION_SHA}" \
  --output "${STAGE}/evaluation" \
  --bootstrap-resamples 10000 \
  --bootstrap-seed 17 \
  > "${STAGE}/validation/evaluation.stdout.json" \
  2> "${STAGE}/validation/evaluation.stderr.txt"
python -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
micromamba deactivate

sha256sum \
  "${PREPARE}" "${BASE_EVALUATOR}" "${EVALUATOR}" \
  tests/unit/test_prepare_sequence_evaluation_matrix.py \
  tests/unit/test_evaluate_sequence_task_native_five_seed_rectangle.py \
  "$0" > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}/src:${ROOT} python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
result = json.loads((root / "evaluation/evaluation.json").read_text(encoding="utf-8"))
if (
    result["status"] != "pass"
    or result["fits"] != 500
    or result["prediction_views"] != 625
    or result["evaluation_groups"] != 25
    or result["partial_rectangle_used"]
    or result["cross_task_family_ranking_performed"]
    or result["test_outcomes_read"]
    or result["sealed_data_read"]
    or result["champion_claim_allowed"]
    or result["universal_claim_allowed"]
):
    raise SystemExit("task-native evaluation terminal contract differs")
print(freeze_tree(root, {
    "artifact_class": "sequence_task_native_five_seed_development_evaluation",
    "dataset_id": "gse296875",
    "status": "pass",
    "fits": 500,
    "prediction_views": 625,
    "evaluation_groups": 25,
    "completion_artifacts_sha256": result["completion_artifacts_sha256"],
    "partial_rectangle_used": False,
    "cross_task_family_ranking_performed": False,
    "test_outcomes_read": False,
    "sealed_data_read": False,
    "champion_claim_allowed": False,
    "universal_claim_allowed": False,
}))
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Task-native five-seed evaluation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
