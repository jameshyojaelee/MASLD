#!/usr/bin/env bash
#SBATCH --job-name=model-data-052
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SOURCE_ARTIFACTS_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
TARGET=${ROOT}/executions/model-data-052-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite study-fold artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${SOURCE_ARTIFACTS_SHA256}" "${SOURCE}/ARTIFACTS.json" \
  | sha256sum --check --strict

module purge
module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
python -m unittest tests.unit.test_build_resource_atlas_study_folds_50000 \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
python scripts/build_resource_atlas_study_folds_50000.py \
  --source "${SOURCE}" \
  --output "${STAGE}/split" \
  --expected-source-artifacts-sha256 "${SOURCE_ARTIFACTS_SHA256}" \
  > "${STAGE}/validation/builder_stdout.json" \
  2> "${STAGE}/validation/builder_stderr.txt"
python - "${STAGE}/split" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

value = verify_frozen_tree(Path(sys.argv[1]))
if (
    value["metadata"].get("rows") != 50_000
    or value["metadata"].get("donors") != 102
    or value["metadata"].get("studies") != 7
    or value["metadata"].get("target_labels_used_for_assignment") is not False
    or value["metadata"].get("sealed_outcomes_read") is not False
):
    raise SystemExit("study-fold frozen receipt differs")
PY
python -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/build_resource_atlas_study_folds_50000.py \
  tests/unit/test_build_resource_atlas_study_folds_50000.py \
  slurm/build_resource_atlas_study_folds_50000.sbatch \
  "${SOURCE}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "cell_state_study_split_campaign",
    "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
    "split_id": "resource_atlas_study_outer_5fold_v1",
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/split/ARTIFACTS.json"
