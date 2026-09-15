#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-diag-505
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
SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SOURCE_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
SPLIT=${ROOT}/executions/model-data-052-21077595/split
SPLIT_SHA256=10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680
TARGET=${ROOT}/executions/model-cpu-diag-505-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite all-fold elastic-net diagnostics.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m py_compile scripts/diagnose_elastic_net_convergence_study_50000.py
"${PYTHON}" -m unittest -v \
  tests.unit.test_diagnose_elastic_net_convergence_study_50000 \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

pids=()
for fold in 0 1 2 3 4; do
  install -d -m 0750 "${STAGE}/fold-${fold}"
  env OMP_NUM_THREADS=3 MKL_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 \
    "${PYTHON}" scripts/diagnose_elastic_net_convergence_study_50000.py \
      --source "${SOURCE}" \
      --source-artifacts-sha256 "${SOURCE_SHA256}" \
      --split "${SPLIT}" \
      --split-artifacts-sha256 "${SPLIT_SHA256}" \
      --outer-fold "${fold}" \
      --config-id float64_tol1e-3_c1_max50000 \
      --config-id float64_tol1e-2_c1 \
      --output "${STAGE}/fold-${fold}/diagnostic" \
      > "${STAGE}/fold-${fold}/stdout.json" \
      2> "${STAGE}/fold-${fold}/stderr.txt" &
  pids+=("$!")
done
for pid in "${pids[@]}"; do
  wait "${pid}"
done

sha256sum \
  scripts/diagnose_elastic_net_convergence_study_50000.py \
  tests/unit/test_diagnose_elastic_net_convergence_study_50000.py \
  slurm/diagnose_elastic_net_all_folds_study_50000.sbatch \
  "${SOURCE}/ARTIFACTS.json" \
  "${SPLIT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
for fold in range(5):
    value = json.loads((root / f"fold-{fold}/diagnostic/diagnostic.json").read_text())
    if (
        value.get("outer_fold") != fold
        or value.get("metrics_calculated") is not False
        or value.get("test_partition_predicted") is not False
        or {item["config_id"] for item in value.get("fits", [])}
        != {"float64_tol1e-3_c1_max50000", "float64_tol1e-2_c1"}
    ):
        raise SystemExit("all-fold diagnostic firewall or roster differs")
freeze_tree(
    root,
    {
        "artifact_class": "cell_baseline_elastic_net_all_fold_convergence_diagnostic",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "outer_folds": 5,
        "logical_tasks": 5,
        "model_selected": False,
        "metrics_calculated": False,
        "test_partition_predicted": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
