#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-train-503
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
MICROMAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
RUNTIME_PREFIX=/gpfs/commons/home/jameslee/micromamba/envs/immune-atlas
RUNTIME_LOCK=${ROOT}/executions/environment-locks/cpu_baseline_smoke--8a7dfd2a15a635d5/environment.lock.json
RUNTIME_LOCK_SHA256=8a7dfd2a15a635d5ef5b25fa8ee68355e6c6cfd047ec464cdc1dea83b7201298
SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SOURCE_ARTIFACTS_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
SPLIT=${ROOT}/executions/model-data-052-21077595/split
SPLIT_ARTIFACTS_SHA256=10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680
DIAGNOSTIC=${ROOT}/executions/model-cpu-diag-507-21084017
DIAGNOSTIC_ARTIFACTS_SHA256=24f2811121a4a31d7e9c77d8fb5aa79e1324c6ef87548d41d2d1714838adba2a
FAILED_STDERR=${ROOT}/executions/model-cpu-train-501-21085652.failed/validation/run_stderr.txt
FAILED_STDERR_SHA256=49834387f973e16e7287bb6c8aaf8f9c821e945b2bab73176134e606b4e51d02
SCRIPT=${ROOT}/scripts/fit_predict_cell_baselines_study_50000.py
SCRIPT_SHA256=66c7d2e9eadac7ca72ef0a7fd9dad92b9f6abb7fd961213a1dff96b03a2bc9df
TEST=${ROOT}/tests/unit/test_fit_predict_cell_baselines_study_50000.py
TEST_SHA256=8bcc299c0905618e77a376ae20c7d1e20df49b6b5c1db2221daad618ea4078de
RESCUE_TEST=${ROOT}/tests/unit/test_sgd_elastic_net_rescue_study_50000.py
RESCUE_TEST_SHA256=2a9291a8f9bbeaae7602b1ae8129d91939caecb8c61d3a42cfeca2bf5d2f90a6
TARGET=${ROOT}/executions/model-cpu-train-503-${SLURM_JOB_ID}
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
  printf 'Refusing to overwrite revised cell-baseline artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
for pair in \
  "${RUNTIME_LOCK_SHA256} ${RUNTIME_LOCK}" \
  "${SOURCE_ARTIFACTS_SHA256} ${SOURCE}/ARTIFACTS.json" \
  "${SPLIT_ARTIFACTS_SHA256} ${SPLIT}/ARTIFACTS.json" \
  "${DIAGNOSTIC_ARTIFACTS_SHA256} ${DIAGNOSTIC}/ARTIFACTS.json" \
  "${FAILED_STDERR_SHA256} ${FAILED_STDERR}" \
  "${SCRIPT_SHA256} ${SCRIPT}" \
  "${TEST_SHA256} ${TEST}" \
  "${RESCUE_TEST_SHA256} ${RESCUE_TEST}"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
if [[ $(find "${ROOT}/executions/model-cpu-train-501-21085652.failed/predictions/folds" \
  -name fold_receipt.json -type f | wc -l) -ne 13 ]]; then
  printf 'Superseded campaign failure census differs.\n' >&2
  exit 1
fi

module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
python scripts/freeze_python_environment.py \
  --micromamba "${MICROMAMBA}" --prefix "${RUNTIME_PREFIX}" \
  --runtime-id cpu_baseline_smoke --output-root "${STAGE}/runtime_probe" \
  > "${STAGE}/validation/runtime_probe.json"
observed_lock=$(find "${STAGE}/runtime_probe" -mindepth 2 -maxdepth 2 \
  -name environment.lock.json -print -quit)
if [[ -z ${observed_lock} ]]; then
  printf 'Runtime probe did not produce an environment lock.\n' >&2
  exit 1
fi
printf '%s  %s\n' "${RUNTIME_LOCK_SHA256}" "${observed_lock}" \
  | sha256sum --check --strict
python - "${SOURCE}" "${SPLIT}" "${DIAGNOSTIC}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
PY
module purge

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
"${MICROMAMBA}" run --prefix "${RUNTIME_PREFIX}" env \
  PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}:${ROOT}/src" \
  OMP_NUM_THREADS="${OMP_NUM_THREADS}" MKL_NUM_THREADS="${MKL_NUM_THREADS}" \
  OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS}" \
  python -m unittest -v \
    tests.unit.test_fit_predict_cell_baselines_study_50000 \
    tests.unit.test_sgd_elastic_net_rescue_study_50000 \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${MICROMAMBA}" run --prefix "${RUNTIME_PREFIX}" env \
  PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}:${ROOT}/src" \
  OMP_NUM_THREADS="${OMP_NUM_THREADS}" MKL_NUM_THREADS="${MKL_NUM_THREADS}" \
  OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS}" \
  python "${SCRIPT}" \
    --source "${SOURCE}" --split "${SPLIT}" \
    --output "${STAGE}/predictions" \
    --expected-source-artifacts-sha256 "${SOURCE_ARTIFACTS_SHA256}" \
    --expected-split-artifacts-sha256 "${SPLIT_ARTIFACTS_SHA256}" \
  > "${STAGE}/validation/run_stdout.json" \
  2> "${STAGE}/validation/run_stderr.txt"
"${MICROMAMBA}" run --prefix "${RUNTIME_PREFIX}" env PYTHONNOUSERSITE=1 \
  python -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/fit_predict_cell_baselines_study_50000.py \
  tests/unit/test_fit_predict_cell_baselines_study_50000.py \
  tests/unit/test_sgd_elastic_net_rescue_study_50000.py \
  slurm/run_cell_baselines_study_50000_v2.sbatch \
  "${RUNTIME_LOCK}" "${SOURCE}/ARTIFACTS.json" \
  "${SPLIT}/ARTIFACTS.json" "${DIAGNOSTIC}/ARTIFACTS.json" \
  "${FAILED_STDERR}" > "${STAGE}/source.sha256"
module load python/3.11.5-GCCcore-13.2.0
python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
predictions = root / "predictions"
receipt = json.loads((predictions / "prediction_receipt.json").read_text())
if (
    receipt.get("status") != "pass_development_predictions"
    or receipt.get("rows") != 50_000
    or receipt.get("donors") != 102
    or receipt.get("studies") != 7
    or len(receipt.get("folds", [])) != 15
    or receipt.get("parameters", {}).get("elastic_net_outer_fit_rule")
    != "inner_loss_rank_until_outer_training_convergence"
    or receipt.get("metrics_calculated") is not False
    or receipt.get("histology_read") is not False
    or receipt.get("sealed_outcomes_read") is not False
):
    raise SystemExit("revised baseline prediction receipt differs")
fallback_ranks = []
for seed in (20260824, 20260825, 20260826):
    for fold in range(5):
        state = json.loads(
            (predictions / f"folds/seed{seed}/fold{fold}/logistic_models.json").read_text()
        )["hvg_pca_elastic_net"]
        if (
            state.get("converged_before_ceiling") is not True
            or state.get("outer_test_partition_used_for_fallback") is not False
            or state.get("selected_inner_validation_rank")
            != len(state.get("outer_training_convergence_fallback", []))
        ):
            raise SystemExit("outer-training convergence fallback receipt differs")
        fallback_ranks.append(int(state["selected_inner_validation_rank"]))
freeze_tree(predictions, {
    "artifact_class": "cell_state_task_native_predictions",
    "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
    "split_id": "resource_atlas_study_outer_5fold_v1",
    "rows": 50_000, "donors": 102, "studies": 7,
    "seeds": [20260824, 20260825, 20260826],
    "outer_training_convergence_rule": "inner_loss_rank_until_convergence",
    "maximum_selected_inner_validation_rank": max(fallback_ranks),
    "prediction_ensemble": "unweighted_mean_probability",
    "metrics_calculated": False, "histology_read": False,
    "sealed_outcomes_read": False, "status": "passed",
})
value = verify_frozen_tree(predictions)
if value["metadata"].get("metrics_calculated") is not False:
    raise SystemExit("baseline prediction artifact differs")
freeze_tree(root, {
    "artifact_class": "cell_state_task_native_campaign",
    "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
    "split_id": "resource_atlas_study_outer_5fold_v1",
    "supersedes_failed_job": "21085652",
    "superseded_failure": "seed20260826_fold3_outer_training_nonconvergence",
    "outer_training_convergence_rule": "inner_loss_rank_until_convergence",
    "maximum_selected_inner_validation_rank": max(fallback_ranks),
    "metrics_calculated": False, "histology_read": False,
    "sealed_outcomes_read": False, "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/predictions/ARTIFACTS.json"
