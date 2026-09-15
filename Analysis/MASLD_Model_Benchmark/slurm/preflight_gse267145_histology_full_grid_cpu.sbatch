#!/usr/bin/env bash
#SBATCH --job-name=model-training-067
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-training-067-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-training-067-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
MOLECULAR=${ROOT}/executions/model-data-064-21079902/fixture/molecular
FOLDS=${ROOT}/executions/model-data-064-21079902/fixture/folds
OUTCOMES=${ROOT}/executions/model-data-061-21079623/activation
SURFACE=${ROOT}/config/evaluation/gse267145_histology_benchmark_surface.json
TARGET=${ROOT}/executions/model-training-067-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite preflight output: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

preserve_failure() {
  local status=$?
  if [[ ${status} -ne 0 && -d "${STAGE}" ]]; then
    printf '%s\n' "${status}" > "${STAGE}/EXIT_STATUS"
    mv "${STAGE}" "${FAILED}"
  fi
  exit "${status}"
}
trap preserve_failure EXIT

cd "${ROOT}"
printf '%s  %s\n' \
  f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb \
  "${RUNTIME}/ARTIFACTS.json" | sha256sum --check --strict
MOLECULAR_SHA=$(sha256sum "${MOLECULAR}/ARTIFACTS.json" | cut -d' ' -f1)
FOLDS_SHA=$(sha256sum "${FOLDS}/ARTIFACTS.json" | cut -d' ' -f1)
OUTCOMES_SHA=$(sha256sum "${OUTCOMES}/ARTIFACTS.json" | cut -d' ' -f1)
SURFACE_SHA=$(sha256sum "${SURFACE}" | cut -d' ' -f1)
[[ "${MOLECULAR_SHA}" == 2a60ef5d5d904f5f833833ebef75e7946e7db3b6b9aaf502cb84ddedf97eb456 ]]
[[ "${FOLDS_SHA}" == 9f5b96e69f217ba643b4b2fd4f3e165049b0cc2bfe9d65c99a7789816b4d7ee3 ]]
[[ "${OUTCOMES_SHA}" == 98a74b97f6a9712a109d2516c3803a4acf15030ed607f861e1ce8c22b4f9753a ]]
sha256sum \
  config/evaluation/gse267145_histology_benchmark_surface.json \
  config/evaluation/gse267145_histology_state_task.toml \
  scripts/build_gse267145_histology_fit_views.py \
  scripts/fit_gse267145_histology_baselines.py \
  scripts/evaluate_gse267145_histology_predictions.py \
  scripts/validate_gse267145_histology_preflight_predictions.py \
  scripts/validate_gse267145_histology_benchmark_surface.py \
  tests/unit/test_validate_gse267145_histology_benchmark_surface.py \
  tests/unit/test_build_gse267145_histology_fit_views.py \
  tests/unit/test_fit_gse267145_histology_baselines.py \
  tests/unit/test_evaluate_gse267145_histology_predictions.py \
  tests/unit/test_validate_gse267145_histology_preflight_predictions.py \
  slurm/preflight_gse267145_histology_full_grid_cpu.sbatch \
  > "${STAGE}/source.sha256"
SOURCE_SHA=$(sha256sum "${STAGE}/source.sha256" | cut -d' ' -f1)

"${PYTHON}" -m py_compile \
  scripts/build_gse267145_histology_fit_views.py \
  scripts/fit_gse267145_histology_baselines.py \
  scripts/evaluate_gse267145_histology_predictions.py \
  scripts/validate_gse267145_histology_preflight_predictions.py \
  scripts/validate_gse267145_histology_benchmark_surface.py
"${PYTHON}" -m unittest -v \
  tests.unit.test_validate_gse267145_histology_benchmark_surface \
  tests.unit.test_build_gse267145_histology_fit_views \
  tests.unit.test_fit_gse267145_histology_baselines \
  tests.unit.test_evaluate_gse267145_histology_predictions \
  tests.unit.test_validate_gse267145_histology_preflight_predictions \
  > "${STAGE}/unit_tests.txt" 2>&1

PYTHONPATH=src "${PYTHON}" scripts/validate_gse267145_histology_benchmark_surface.py \
  --root "${ROOT}" \
  --surface "${SURFACE}" \
  --output "${STAGE}/design_validation"
PYTHONPATH=src "${PYTHON}" - "${STAGE}/design_validation" "${SURFACE_SHA}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
print(freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "gse267145_histology_updated_surface_validation",
    "surface_sha256": sys.argv[2],
    "status": "passed_design_only_no_fit",
}))
PY
DESIGN_SHA=$(sha256sum "${STAGE}/design_validation/ARTIFACTS.json" | cut -d' ' -f1)

"${PYTHON}" scripts/build_gse267145_histology_fit_views.py \
  --molecular "${MOLECULAR}" \
  --molecular-artifacts-sha256 "${MOLECULAR_SHA}" \
  --folds "${FOLDS}" \
  --folds-artifacts-sha256 "${FOLDS_SHA}" \
  --outcomes "${OUTCOMES}" \
  --outcomes-artifacts-sha256 "${OUTCOMES_SHA}" \
  --inner-validation "${STAGE}/design_validation" \
  --inner-validation-artifacts-sha256 "${DESIGN_SHA}" \
  --output "${STAGE}/fit_views"
PYTHONPATH=src "${PYTHON}" - "${STAGE}/fit_views" "${DESIGN_SHA}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
print(freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "gse267145_histology_fold_scoped_fit_views",
    "inner_validation_artifacts_sha256": sys.argv[2],
    "held_outcomes_included": False,
    "status": "passed",
}))
PY
FIT_VIEWS_SHA=$(sha256sum "${STAGE}/fit_views/ARTIFACTS.json" | cut -d' ' -f1)

/usr/bin/time \
  -f '{"elapsed_seconds": %e, "maximum_resident_kib": %M, "user_seconds": %U, "system_seconds": %S}' \
  -o "${STAGE}/fit_resource.json" \
  "${PYTHON}" scripts/fit_gse267145_histology_baselines.py \
    --molecular "${MOLECULAR}" \
    --molecular-artifacts-sha256 "${MOLECULAR_SHA}" \
    --folds "${FOLDS}" \
    --folds-artifacts-sha256 "${FOLDS_SHA}" \
    --fit-views "${STAGE}/fit_views" \
    --fit-views-artifacts-sha256 "${FIT_VIEWS_SHA}" \
    --surface "${SURFACE}" \
    --surface-sha256 "${SURFACE_SHA}" \
    --preflight-outer-fold 0 \
    --preflight-seed 1701 \
    --output "${STAGE}/preflight_fit" \
    > "${STAGE}/fit_stdout.json"

PYTHONPATH=src "${PYTHON}" - "${STAGE}/preflight_fit/predictions" "${SURFACE_SHA}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
print(freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "gse267145_histology_full_grid_timing_predictions",
    "surface_sha256": sys.argv[2],
    "outer_fold": 0,
    "seed": 1701,
    "prediction_files": 11,
    "outcomes_read": False,
    "metrics_calculated": False,
    "status": "passed_unscored",
}))
PY
PREDICTIONS_SHA=$(sha256sum "${STAGE}/preflight_fit/predictions/ARTIFACTS.json" | cut -d' ' -f1)

"${PYTHON}" scripts/validate_gse267145_histology_preflight_predictions.py \
  --predictions "${STAGE}/preflight_fit/predictions" \
  --predictions-artifacts-sha256 "${PREDICTIONS_SHA}" \
  --outcomes "${OUTCOMES}" \
  --outcomes-artifacts-sha256 "${OUTCOMES_SHA}" \
  --folds "${FOLDS}" \
  --folds-artifacts-sha256 "${FOLDS_SHA}" \
  --outer-fold 0 \
  --seed 1701 \
  --output "${STAGE}/evaluator_validation" \
  > "${STAGE}/evaluator_stdout.json"
PYTHONPATH=src "${PYTHON}" - "${STAGE}/evaluator_validation" "${PREDICTIONS_SHA}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
print(freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "gse267145_histology_preflight_evaluator_validation",
    "predictions_artifacts_sha256": sys.argv[2],
    "metrics_calculated": False,
    "status": "passed_firewall_validation_no_metrics",
}))
PY

"${PYTHON}" -m pip freeze > "${STAGE}/pip_freeze.txt"
PYTHONPATH=src "${PYTHON}" - \
  "${STAGE}" "${SURFACE_SHA}" "${DESIGN_SHA}" "${FIT_VIEWS_SHA}" "${PREDICTIONS_SHA}" "${SOURCE_SHA}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, write_json_exclusive
root = Path(sys.argv[1])
fit = json.loads((root / "preflight_fit/receipt.json").read_text())
resource = json.loads((root / "fit_resource.json").read_text())
write_json_exclusive(root / "preflight_receipt.json", {
    "schema_version": "masld-bench-gse267145-full-grid-timing-preflight-v1",
    "status": "passed",
    "surface_sha256": sys.argv[2],
    "design_validation_artifacts_sha256": sys.argv[3],
    "fit_views_artifacts_sha256": sys.argv[4],
    "predictions_artifacts_sha256": sys.argv[5],
    "source_manifest_sha256": sys.argv[6],
    "outer_fold": 0,
    "seed": 1701,
    "logical_outer_seed_tasks": 1,
    "model_family_fits": 11,
    "exact_production_grids_used": True,
    "fit_counts": fit["fit_counts"],
    "invalid_candidate_count": fit["invalid_candidate_count"],
    "invalid_candidate_reasons": fit["invalid_candidate_reasons"],
    "resource_measurement": resource,
    "production_predictions_complete": False,
    "metrics_calculated": False,
})
print(freeze_tree(root, {
    "artifact_class": "gse267145_histology_full_grid_timing_preflight",
    "surface_sha256": sys.argv[2],
    "source_manifest_sha256": sys.argv[6],
    "outer_fold": 0,
    "seed": 1701,
    "model_family_fits": 11,
    "metrics_calculated": False,
    "status": "passed",
}))
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'GSE267145 full-grid timing preflight: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
