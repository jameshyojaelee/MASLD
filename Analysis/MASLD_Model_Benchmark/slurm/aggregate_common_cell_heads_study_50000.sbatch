#!/usr/bin/env bash
#SBATCH --job-name=model-aggregate
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
: "${SHARD_MANIFEST:?absolute shard manifest is required}"
: "${SHARD_MANIFEST_SHA256:?shard manifest SHA-256 is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SOURCE_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
SPLIT=${ROOT}/executions/model-data-052-21077595/split
SPLIT_SHA256=10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680
TARGET=${ROOT}/executions/model-aggregate-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite common-head aggregate.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SOURCE_SHA256}" "${SOURCE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SPLIT_SHA256}" "${SPLIT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SHARD_MANIFEST_SHA256}" "${SHARD_MANIFEST}" \
  | sha256sum --check --strict

cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
"${PYTHON}" -m unittest -v \
  tests.unit.test_fit_predict_common_cell_heads_study_50000 \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${PYTHON}" scripts/fit_predict_common_cell_heads_study_50000.py aggregate \
  --shard-manifest "${SHARD_MANIFEST}" \
  --source "${SOURCE}" \
  --split "${SPLIT}" \
  --output "${STAGE}/aggregate" \
  --expected-source-artifacts-sha256 "${SOURCE_SHA256}" \
  --expected-split-artifacts-sha256 "${SPLIT_SHA256}" \
  > "${STAGE}/validation/aggregate.stdout.json" \
  2> "${STAGE}/validation/aggregate.stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/fit_predict_common_cell_heads_study_50000.py \
  tests/unit/test_fit_predict_common_cell_heads_study_50000.py \
  slurm/aggregate_common_cell_heads_study_50000.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  "${SOURCE}/ARTIFACTS.json" \
  "${SPLIT}/ARTIFACTS.json" \
  "${SHARD_MANIFEST}" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
value = verify_frozen_tree(root / "aggregate")
receipt = json.loads((root / "aggregate/prediction_receipt.json").read_text())
if (
    value["metadata"].get("logical_shards") != 30
    or receipt.get("logical_shards") != 30
    or receipt.get("checkpoint_selection_uses_validation_labels") is not False
    or receipt.get("inner_models_used_for_outer_prediction") is not False
    or receipt.get("outer_prediction_uses_final_refit") is not True
    or receipt.get("prediction_tables_contain_observed_labels") is not False
    or receipt.get("sealed_outcomes_read") is not False
):
    raise SystemExit("common-head aggregate receipt differs")
freeze_tree(root, {
    "artifact_class": "cell_foundation_common_head_aggregate_campaign",
    "logical_shards": 30,
    "metrics_calculated": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/aggregate/ARTIFACTS.json"
