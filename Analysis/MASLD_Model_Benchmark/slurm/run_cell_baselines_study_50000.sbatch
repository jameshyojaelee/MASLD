#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-train-501
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark
MICROMAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
RUNTIME_PREFIX=/gpfs/commons/home/jameslee/micromamba/envs/immune-atlas
RUNTIME_LOCK=${ROOT}/executions/environment-locks/cpu_baseline_smoke--8a7dfd2a15a635d5/environment.lock.json
RUNTIME_LOCK_SHA256=8a7dfd2a15a635d5ef5b25fa8ee68355e6c6cfd047ec464cdc1dea83b7201298
SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SOURCE_ARTIFACTS_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
SPLIT=${ROOT}/executions/model-data-052-21077595/split
SPLIT_ARTIFACTS_SHA256=10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680
TARGET=${ROOT}/executions/model-cpu-train-501-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite cell-baseline artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_LOCK_SHA256}" "${RUNTIME_LOCK}" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SOURCE_ARTIFACTS_SHA256}" "${SOURCE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SPLIT_ARTIFACTS_SHA256}" "${SPLIT}/ARTIFACTS.json" \
  | sha256sum --check --strict

module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
cd "${ROOT}"
python scripts/freeze_python_environment.py \
  --micromamba "${MICROMAMBA}" \
  --prefix "${RUNTIME_PREFIX}" \
  --runtime-id cpu_baseline_smoke \
  --output-root "${STAGE}/runtime_probe" \
  > "${STAGE}/validation/runtime_probe.json"
observed_lock=$(find "${STAGE}/runtime_probe" -mindepth 2 -maxdepth 2 -name environment.lock.json -print -quit)
if [[ -z "${observed_lock}" ]]; then
  printf 'Runtime probe did not produce an environment lock.\n' >&2
  exit 1
fi
printf '%s  %s\n' "${RUNTIME_LOCK_SHA256}" "${observed_lock}" \
  | sha256sum --check --strict
python - "${SOURCE}" "${SPLIT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
verify_frozen_tree(Path(sys.argv[2]))
PY
module purge

export PYTHONPATH="${ROOT}:${ROOT}/src"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
"${MICROMAMBA}" run --prefix "${RUNTIME_PREFIX}" env \
  PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}:${ROOT}/src" \
  OMP_NUM_THREADS="${OMP_NUM_THREADS}" MKL_NUM_THREADS="${MKL_NUM_THREADS}" \
  OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS}" \
  python -m unittest tests.unit.test_fit_predict_cell_baselines_study_50000 \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${MICROMAMBA}" run --prefix "${RUNTIME_PREFIX}" env \
  PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}:${ROOT}/src" \
  OMP_NUM_THREADS="${OMP_NUM_THREADS}" MKL_NUM_THREADS="${MKL_NUM_THREADS}" \
  OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS}" \
  python scripts/fit_predict_cell_baselines_study_50000.py \
  --source "${SOURCE}" \
  --split "${SPLIT}" \
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
  slurm/run_cell_baselines_study_50000.sbatch \
  "${RUNTIME_LOCK}" \
  "${SOURCE}/ARTIFACTS.json" \
  "${SPLIT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
module load python/3.11.5-GCCcore-13.2.0
python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1])
freeze_tree(root / "predictions", {
    "artifact_class": "cell_state_task_native_predictions",
    "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
    "split_id": "resource_atlas_study_outer_5fold_v1",
    "rows": 50_000,
    "donors": 102,
    "studies": 7,
    "seeds": [20260824, 20260825, 20260826],
    "prediction_ensemble": "unweighted_mean_probability",
    "metrics_calculated": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
value = verify_frozen_tree(root / "predictions")
if (
    value["metadata"].get("rows") != 50_000
    or value["metadata"].get("donors") != 102
    or value["metadata"].get("studies") != 7
    or value["metadata"].get("metrics_calculated") is not False
    or value["metadata"].get("sealed_outcomes_read") is not False
):
    raise SystemExit("baseline prediction receipt differs")
freeze_tree(root, {
    "artifact_class": "cell_state_task_native_campaign",
    "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
    "split_id": "resource_atlas_study_outer_5fold_v1",
    "metrics_calculated": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/predictions/ARTIFACTS.json"
