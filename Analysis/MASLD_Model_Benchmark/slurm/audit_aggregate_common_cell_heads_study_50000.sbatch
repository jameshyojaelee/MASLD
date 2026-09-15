#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-audit-503
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${BUNDLE_ROOT:?BUNDLE_ROOT is required}"
: "${BUNDLE_SHA256:?BUNDLE_SHA256 is required}"
: "${EXPECTED_MODEL_ID:?EXPECTED_MODEL_ID is required}"

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
SOURCE_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
SPLIT=${ROOT}/executions/model-data-052-21077595/split
SPLIT_SHA256=10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680
TARGET=${ROOT}/executions/model-cpu-audit-503-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite 50,000-cell audit artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${BUNDLE_SHA256}" "${BUNDLE_ROOT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SPLIT_SHA256}" "${SPLIT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict

cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
"${PYTHON}" -m unittest -v \
  tests.unit.test_audit_evaluate_common_cell_heads_study_50000 \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${PYTHON}" scripts/audit_aggregate_common_cell_heads_study_50000.py \
  --bundle "${BUNDLE_ROOT}" \
  --split "${SPLIT}" \
  --output "${STAGE}/predictions" \
  --expected-bundle-artifacts-sha256 "${BUNDLE_SHA256}" \
  --expected-split-artifacts-sha256 "${SPLIT_SHA256}" \
  --expected-source-artifacts-sha256 "${SOURCE_SHA256}" \
  --expected-model-id "${EXPECTED_MODEL_ID}" \
  > "${STAGE}/validation/aggregator_stdout.json" \
  2> "${STAGE}/validation/aggregator_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/audit_aggregate_common_cell_heads_study_50000.py \
  scripts/score_common_cell_heads_study_50000.py \
  scripts/score_cell_baselines_study_50000.py \
  tests/unit/test_audit_evaluate_common_cell_heads_study_50000.py \
  slurm/audit_aggregate_common_cell_heads_study_50000.sbatch \
  slurm/score_common_cell_heads_study_50000.sbatch \
  "${BUNDLE_ROOT}/ARTIFACTS.json" \
  "${SPLIT}/ARTIFACTS.json" \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" "${EXPECTED_MODEL_ID}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
expected_model_id = sys.argv[2]
verify_frozen_tree(root / "predictions")
receipt = json.loads((root / "predictions/prediction_receipt.json").read_text())
if (
    receipt["status"] != "pass_development_predictions_audited"
    or receipt["model_id"] != expected_model_id
    or receipt["rows"] != 50_000
    or receipt["logical_shards"] != 30
    or receipt["prediction_tables"] != 6
    or not receipt["donor_safe_outer_folds_verified"]
    or not receipt["whole_study_outer_folds_verified"]
    or not receipt["head_seed_fold_census_complete"]
    or receipt["metrics_calculated"]
    or receipt["development_labels_read_by_aggregator"]
    or receipt["sealed_outcomes_read"]
):
    raise SystemExit("50,000-cell audited prediction receipt differs")
freeze_tree(root, {
    "artifact_class": "cell_foundation_common_head_audit_campaign",
    "model_id": expected_model_id,
    "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
    "split_id": "resource_atlas_study_outer_5fold_v1",
    "logical_shards": 30,
    "metrics_calculated": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/predictions/ARTIFACTS.json"
