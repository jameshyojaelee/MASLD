#!/usr/bin/env bash
#SBATCH --job-name=model-check
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-check-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-check-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET=${ROOT}/executions/sequence-regulatory-completion-gate-tests-${SLURM_JOB_ID}
STAGE=${TARGET}.staging

preserve_failed_stage() {
  local failed=${TARGET}.failed
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite validation target: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
printf '%s  %s\n' \
  edede235b0b87913de334142a25d184a89803870cb92c640a54a81375d68ec08 config/evaluation/sequence_regulatory_full_rectangle_completion_20260825.json \
  69a7f5a504476482b2e32a9486a9f37b43b4e3cb3cf297422026303be03b39da config/evaluation/scbasset_five_seed_rectangular_evaluator_20260825.json \
  afba2e8882c84596451a0d74ac387b6ca587bdbde825c35352b5246ade8ca491 scripts/audit_sequence_regulatory_full_rectangle_completion.py \
  4ce7a88a08c1da6c316b85d8a15b32e637cea4243d0ffadac8883fef79195b12 scripts/evaluate_sequence_task_native_five_seed_rectangle.py \
  1c314b6aee1ec0d8043fd8a1221c8ee69b5781bb240acbf78d4b8a30feb353e6 scripts/evaluate_scbasset_five_seed_rectangle.py \
  4d69663edc21281eacad6e3ce8b5c590f696a70261da07cd71fe6b5eb95c65e6 tests/unit/test_sequence_regulatory_full_rectangle_completion.py \
  0f04f22378eed5a794e7546889bd5d1d2443fb36a3294aac8a6db69e4d3ec589 tests/unit/test_evaluate_sequence_task_native_five_seed_rectangle.py \
  7c230d85f1caafff06a2fc1be7bf79a5e5dc3857b2f93d1833f17893e8f37027 tests/unit/test_evaluate_scbasset_five_seed_rectangle.py \
  | sha256sum --check --strict

PYTHONPATH=${ROOT}/src:${ROOT} python -m unittest -v \
  tests.unit.test_sequence_regulatory_full_rectangle_completion \
  tests.unit.test_evaluate_sequence_task_native_five_seed_rectangle \
  tests.unit.test_evaluate_scbasset_five_seed_rectangle \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
python --version > "${STAGE}/validation/python-version.txt"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  config/evaluation/sequence_regulatory_full_rectangle_completion_20260825.json \
  config/evaluation/scbasset_five_seed_rectangular_evaluator_20260825.json \
  scripts/audit_sequence_regulatory_full_rectangle_completion.py \
  scripts/evaluate_sequence_task_native_five_seed_rectangle.py \
  scripts/evaluate_scbasset_five_seed_rectangle.py \
  tests/unit/test_sequence_regulatory_full_rectangle_completion.py \
  tests/unit/test_evaluate_sequence_task_native_five_seed_rectangle.py \
  tests/unit/test_evaluate_scbasset_five_seed_rectangle.py \
  "$0" > "${STAGE}/source.sha256"

PYTHONPATH=${ROOT}/src:${ROOT} python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
stdout = (root / "validation/tests.stdout.txt").read_text(encoding="utf-8")
stderr = (root / "validation/tests.stderr.txt").read_text(encoding="utf-8")
if "OK" not in stderr or "FAILED" in stderr or "ERROR" in stderr:
    raise SystemExit("completion-gate unit tests did not pass")
print(freeze_tree(root, {
    "artifact_class": "sequence_regulatory_completion_gate_acceptance",
    "status": "pass",
    "unit_tests_passed": True,
    "outcomes_read": False,
    "prediction_values_read": False,
    "sealed_data_read": False,
    "partial_ranking_performed": False,
    "source_stdout_nonempty": bool(stdout.strip()),
}))
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Completion-gate acceptance: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
