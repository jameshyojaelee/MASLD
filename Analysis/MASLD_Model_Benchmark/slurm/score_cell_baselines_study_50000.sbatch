#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-score-502
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
: "${PREDICTIONS_ROOT:?PREDICTIONS_ROOT is required}"
: "${PREDICTIONS_SHA256:?PREDICTIONS_SHA256 is required}"

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SOURCE_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
SPLIT=${ROOT}/executions/model-data-052-21077595/split
SPLIT_SHA256=10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680
TARGET=${ROOT}/executions/model-cpu-score-502-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite 50,000-cell scorer artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${SOURCE_SHA256}" "${SOURCE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SPLIT_SHA256}" "${SPLIT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${PREDICTIONS_SHA256}" "${PREDICTIONS_ROOT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict

cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
"${PYTHON}" -m unittest -v tests.unit.test_score_cell_baselines_study_50000 \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${PYTHON}" scripts/score_cell_baselines_study_50000.py \
  --source "${SOURCE}" \
  --split "${SPLIT}" \
  --predictions "${PREDICTIONS_ROOT}" \
  --output "${STAGE}/evaluation" \
  --expected-source-artifacts-sha256 "${SOURCE_SHA256}" \
  --expected-split-artifacts-sha256 "${SPLIT_SHA256}" \
  --expected-prediction-artifacts-sha256 "${PREDICTIONS_SHA256}" \
  > "${STAGE}/validation/evaluator_stdout.json" \
  2> "${STAGE}/validation/evaluator_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/score_cell_baselines_study_50000.py \
  tests/unit/test_score_cell_baselines_study_50000.py \
  slurm/score_cell_baselines_study_50000.sbatch \
  "${SOURCE}/ARTIFACTS.json" \
  "${SPLIT}/ARTIFACTS.json" \
  "${PREDICTIONS_ROOT}/ARTIFACTS.json" \
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
    or receipt["rows"] != 50_000
    or receipt["donors"] != 102
    or receipt["studies"] != 7
    or receipt["models_fit"]
    or receipt["preprocessing_fit"]
    or receipt["calibration_fit"]
    or receipt["histology_read"]
    or receipt["sealed_outcomes_read"]
    or receipt["bootstrap"]["resamples"] != 10_000
    or receipt["bootstrap"]["cells_are_independent_units"]
):
    raise SystemExit("50,000-cell evaluation-only receipt differs")
freeze_tree(root, {
    "artifact_class": "cell_state_study_heldout_evaluation_campaign",
    "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
    "split_id": "resource_atlas_study_outer_5fold_v1",
    "models_fit": False,
    "histology_read": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/evaluation/ARTIFACTS.json"
