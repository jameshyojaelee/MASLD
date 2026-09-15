#!/usr/bin/env bash
#SBATCH --job-name=model-preflight-146
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/scvi-tools-1.4.2-inductive-portable-21079606
RUNTIME_ARTIFACTS_SHA256=d885ada24fdad9f6e3d57a92ff870ca8ff05f07701ac6b0762f55a8cb19db61a
SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SOURCE_ARTIFACTS_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
SPLIT=${ROOT}/executions/model-data-052-21077595/split
SPLIT_ARTIFACTS_SHA256=10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680
WRAPPER=${ROOT}/slurm/run_scvi_scanvi_study_50000_seed_bundle.sbatch
WRAPPER_SHA256=c4661c4d4313419bad9e96f1b06bdbeae3779e7600c0b18c8422f2b4620f5dbe
TARGET=${ROOT}/executions/model-preflight-146-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite scVI/scANVI preflight.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SOURCE_ARTIFACTS_SHA256}" "${SOURCE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SPLIT_ARTIFACTS_SHA256}" "${SPLIT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${WRAPPER_SHA256}" "${WRAPPER}" \
  | sha256sum --check --strict

module purge
module load python/3.12.3-GCCcore-13.3.0
cd "${ROOT}"
export PYTHONPATH=${RUNTIME}/site-packages:${ROOT}:${ROOT}/src
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
python -m unittest \
  tests.unit.test_scvi_scanvi_inductive \
  tests.unit.test_fit_predict_scvi_scanvi_study_50000 \
  tests.unit.test_aggregate_scvi_scanvi_study_50000 \
  > "${STAGE}/validation/unit-tests.txt" \
  2>&1
python scripts/preflight_scvi_scanvi_study_50000.py \
  --source "${SOURCE}" \
  --split "${SPLIT}" \
  --runtime "${RUNTIME}" \
  --wrapper "${WRAPPER}" \
  --output "${STAGE}/preflight" \
  > "${STAGE}/validation/preflight.stdout.json" \
  2> "${STAGE}/validation/preflight.stderr.txt"
cp "${RUNTIME}/validation/pip-freeze.txt" "${STAGE}/validation/runtime-pip-freeze.txt"
sha256sum \
  scripts/preflight_scvi_scanvi_study_50000.py \
  scripts/fit_predict_scvi_scanvi_study_50000.py \
  scripts/aggregate_scvi_scanvi_study_50000.py \
  scripts/score_scvi_scanvi_study_50000.py \
  src/masld_bench/adapters/scvi_scanvi_inductive.py \
  tests/unit/test_scvi_scanvi_inductive.py \
  tests/unit/test_fit_predict_scvi_scanvi_study_50000.py \
  tests/unit/test_aggregate_scvi_scanvi_study_50000.py \
  slurm/preflight_scvi_scanvi_study_50000.sbatch \
  slurm/score_scvi_scanvi_study_50000.sbatch \
  "${WRAPPER}" \
  "${RUNTIME}/ARTIFACTS.json" \
  "${SOURCE}/ARTIFACTS.json" \
  "${SPLIT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}/src python - "${STAGE}" "${RUNTIME_ARTIFACTS_SHA256}" "${WRAPPER_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "preflight/preflight.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "pass_build_only_gpu_pending"
    or receipt.get("rows") != 50_000
    or receipt.get("donors") != 102
    or receipt.get("studies") != 7
    or receipt.get("runtime_artifacts_sha256") != sys.argv[2]
    or receipt.get("gpu_probe_pending") is not True
    or receipt.get("query_training_allowed") is not False
    or receipt.get("query_adaptation_allowed") is not False
    or receipt.get("metrics_calculated") is not False
    or receipt.get("sealed_outcomes_read") is not False
):
    raise SystemExit("scVI/scANVI preflight receipt differs")
freeze_tree(
    root / "preflight",
    {
        "artifact_class": "scvi_scanvi_50000_build_preflight",
        "runtime_artifacts_sha256": sys.argv[2],
        "gpu_wrapper_sha256": sys.argv[3],
        "gpu_probe_pending": True,
        "project_data_model_fit_executed": False,
        "metrics_calculated": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
verify_frozen_tree(root / "preflight")
freeze_tree(
    root,
    {
        "artifact_class": "scvi_scanvi_50000_build_preflight_campaign",
        "runtime_artifacts_sha256": sys.argv[2],
        "gpu_wrapper_sha256": sys.argv[3],
        "gpu_probe_pending": True,
        "project_data_model_fit_executed": False,
        "metrics_calculated": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/preflight/ARTIFACTS.json"
