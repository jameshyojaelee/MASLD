#!/usr/bin/env bash
#SBATCH --job-name=model-data-066
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
MOLECULAR=${ROOT}/executions/model-data-064-21079902/fixture/molecular
FOLDS=${ROOT}/executions/model-data-064-21079902/fixture/folds
OUTCOMES=${ROOT}/executions/model-data-061-21079623/activation
SURFACE=${ROOT}/config/evaluation/gse267145_histology_benchmark_surface.json
TARGET=${ROOT}/executions/model-data-066-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite histology benchmark-surface validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}/src:${ROOT}"

python -m unittest -v tests.unit.test_validate_gse267145_histology_benchmark_surface \
  > "${STAGE}/unit_tests.txt" 2>&1
python - "${MOLECULAR}" "${FOLDS}" "${OUTCOMES}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
print("input frozen-tree verification passed")
PY
python scripts/validate_gse267145_histology_benchmark_surface.py \
  --root "${ROOT}" \
  --surface "${SURFACE}" \
  --output "${STAGE}/validation" \
  > "${STAGE}/run_stdout.json" 2> "${STAGE}/run_stderr.txt"
python -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  "${SURFACE}" \
  scripts/validate_gse267145_histology_benchmark_surface.py \
  tests/unit/test_validate_gse267145_histology_benchmark_surface.py \
  slurm/validate_gse267145_histology_benchmark_surface_cpu.sbatch \
  "${MOLECULAR}/ARTIFACTS.json" \
  "${FOLDS}/ARTIFACTS.json" \
  "${OUTCOMES}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/validation_receipt.json").read_text())
census = json.loads((root / "validation/inner_fold_endpoint_census.json").read_text())
if (
    receipt.get("status") != "passed_design_only_no_fit"
    or receipt.get("baseline_count") != 7
    or receipt.get("production_fit_authorized") is not False
    or receipt.get("model_fitted") is not False
    or receipt.get("metrics_calculated") is not False
    or receipt.get("one_standard_error_selection") is not True
    or receipt.get("endpoint_safe_inner_splits") is not True
    or census.get("every_inner_training_endpoint_safe") is not True
    or census.get("minimum_inner_training_fibrosis3", 0) < 1
    or census.get("minimum_inner_training_threshold_class_count", 0) < 1
    or census.get("algorithm_does_not_use") != ["source_stage5", "recorded_sex"]
):
    raise SystemExit("histology benchmark-surface receipt differs")
freeze_tree(
    root / "validation",
    {
        "artifact_class": "gse267145_histology_benchmark_surface_validation",
        "participants": 99,
        "lanes": 3,
        "baseline_count": 7,
        "production_fit_authorized": False,
        "model_fitted": False,
        "metrics_calculated": False,
        "one_standard_error_selection": True,
        "endpoint_safe_inner_splits": True,
        "source_stage5_selection": False,
        "status": "passed_design_only_no_fit",
    },
)
verify_frozen_tree(root / "validation")
freeze_tree(
    root,
    {
        "artifact_class": "gse267145_histology_benchmark_surface_campaign",
        "production_fit_authorized": False,
        "model_fitted": False,
        "metrics_calculated": False,
        "status": "passed_design_only_no_fit",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/validation/ARTIFACTS.json"
