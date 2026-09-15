#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-diag-507
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/immune-atlas/bin/python
RUNTIME_LOCK=${ROOT}/executions/environment-locks/cpu_baseline_smoke--8a7dfd2a15a635d5/environment.lock.json
RUNTIME_LOCK_SHA256=8a7dfd2a15a635d5ef5b25fa8ee68355e6c6cfd047ec464cdc1dea83b7201298
SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SOURCE_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
SPLIT=${ROOT}/executions/model-data-052-21077595/split
SPLIT_SHA256=10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680
DIAGNOSTIC_504=${ROOT}/executions/model-cpu-diag-504-21079655
DIAGNOSTIC_504_SHA256=bbccba9f03aaa9736a9ac9146e25d00eb7c9fa9b6b2a06c11ce62e4c6b26306e
DIAGNOSTIC_505=${ROOT}/executions/model-cpu-diag-505-21080294
DIAGNOSTIC_505_SHA256=8f00036292a60a1712e4da0f5dd8f9e9c4e8a47541e4af8a15778d145d2124c0
TARGET=${ROOT}/executions/model-cpu-diag-507-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  status=$?
  if [[ ${status} -ne 0 && -d ${STAGE} && ! -e ${FAILED} ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
  exit "${status}"
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite SGD elastic-net rescue diagnostics.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
printf '%s  %s\n' "${RUNTIME_LOCK_SHA256}" "${RUNTIME_LOCK}" | sha256sum --check --strict
printf '%s  %s\n' "${SOURCE_SHA256}" "${SOURCE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${SPLIT_SHA256}" "${SPLIT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${DIAGNOSTIC_504_SHA256}" "${DIAGNOSTIC_504}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${DIAGNOSTIC_505_SHA256}" "${DIAGNOSTIC_505}/ARTIFACTS.json" | sha256sum --check --strict

"${PYTHON}" -m py_compile \
  scripts/fit_predict_cell_baselines_study_50000.py \
  scripts/diagnose_sgd_elastic_net_rescue_study_50000.py
"${PYTHON}" -m unittest -v \
  tests.unit.test_fit_predict_cell_baselines_study_50000 \
  tests.unit.test_sgd_elastic_net_rescue_study_50000 \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

pids=()
for fold in 0 1 2 3 4; do
  install -d -m 0750 "${STAGE}/fold-${fold}"
  env OMP_NUM_THREADS=3 MKL_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 \
    "${PYTHON}" scripts/diagnose_sgd_elastic_net_rescue_study_50000.py \
      --source "${SOURCE}" \
      --source-artifacts-sha256 "${SOURCE_SHA256}" \
      --split "${SPLIT}" \
      --split-artifacts-sha256 "${SPLIT_SHA256}" \
      --diagnostic-504 "${DIAGNOSTIC_504}" \
      --diagnostic-505 "${DIAGNOSTIC_505}" \
      --outer-fold "${fold}" \
      --seed 20260824 \
      --output "${STAGE}/fold-${fold}/diagnostic" \
      > "${STAGE}/fold-${fold}/stdout.json" \
      2> "${STAGE}/fold-${fold}/stderr.txt" &
  pids+=("$!")
done
for pid in "${pids[@]}"; do
  wait "${pid}"
done
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/fit_predict_cell_baselines_study_50000.py \
  scripts/diagnose_sgd_elastic_net_rescue_study_50000.py \
  tests/unit/test_fit_predict_cell_baselines_study_50000.py \
  tests/unit/test_sgd_elastic_net_rescue_study_50000.py \
  slurm/diagnose_sgd_elastic_net_rescue_study_50000.sbatch \
  "${RUNTIME_LOCK}" \
  "${SOURCE}/ARTIFACTS.json" "${SPLIT}/ARTIFACTS.json" \
  "${DIAGNOSTIC_504}/ARTIFACTS.json" "${DIAGNOSTIC_505}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
for fold in range(5):
    value = json.loads(
        (root / f"fold-{fold}/diagnostic/diagnostic.json").read_text(encoding="utf-8")
    )
    if (
        value.get("status") != "passed_training_only_rescue"
        or value.get("outer_fold") != fold
        or value.get("final_fit", {}).get("converged_before_ceiling") is not True
        or value.get("outer_test_features_transformed") is not False
        or value.get("outer_test_partition_predicted") is not False
        or value.get("outer_test_metrics_calculated") is not False
        or value.get("sealed_outcomes_read") is not False
    ):
        raise SystemExit(f"SGD elastic-net rescue diagnostic differs for fold {fold}")
freeze_tree(
    root,
    {
        "artifact_class": "cell_baseline_sgd_elastic_net_all_fold_rescue_diagnostic",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "outer_folds": 5,
        "logical_tasks": 5,
        "metrics_calculated_on_outer_test": False,
        "outer_test_partition_predicted": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
verify_frozen_tree(root)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
