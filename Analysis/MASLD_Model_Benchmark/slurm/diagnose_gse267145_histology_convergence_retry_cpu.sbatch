#!/usr/bin/env bash
#SBATCH --job-name=model-check-070
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=3
export MKL_NUM_THREADS=3
export OPENBLAS_NUM_THREADS=3

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
MOLECULAR=${ROOT}/executions/model-data-064-21079902/fixture/molecular
FOLDS=${ROOT}/executions/model-data-064-21079902/fixture/folds
FIT_VIEWS=${ROOT}/executions/model-training-067-21080776/fit_views
SURFACE=${ROOT}/config/evaluation/gse267145_histology_benchmark_surface.json
TARGET=${ROOT}/executions/model-check-070-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE267145 convergence diagnostics.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${STAGE}/units"
cd "${ROOT}"
export PYTHONPATH="${ROOT}/src:${ROOT}"
"${PYTHON}" -m py_compile scripts/fit_gse267145_histology_baselines.py
"${PYTHON}" -m unittest -v tests.unit.test_fit_gse267145_histology_baselines \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

declare -a units=("1 1701")
declare -a pids=()
for unit in "${units[@]}"; do
  read -r outer seed <<< "${unit}"
  unit_root=${STAGE}/units/outer_${outer}--seed_${seed}
  install -d -m 0750 "${unit_root}"
  /usr/bin/time -f \
    '{"elapsed_seconds": %e, "maximum_resident_kib": %M, "user_seconds": %U, "system_seconds": %S}' \
    -o "${unit_root}/resource.json" \
    "${PYTHON}" scripts/fit_gse267145_histology_baselines.py \
      --molecular "${MOLECULAR}" \
      --molecular-artifacts-sha256 2a60ef5d5d904f5f833833ebef75e7946e7db3b6b9aaf502cb84ddedf97eb456 \
      --folds "${FOLDS}" \
      --folds-artifacts-sha256 9f5b96e69f217ba643b4b2fd4f3e165049b0cc2bfe9d65c99a7789816b4d7ee3 \
      --fit-views "${FIT_VIEWS}" \
      --fit-views-artifacts-sha256 8254274f1b11089656f4b6d5ec446a6ba491e2f28ae8abcc088fc6e78ac75d6a \
      --surface "${SURFACE}" \
      --surface-sha256 23367a28bd16266c79adde854534b124cfad7d2062788625839e78f3267a8c0b \
      --preflight-outer-fold "${outer}" \
      --preflight-seed "${seed}" \
      --output "${unit_root}/fit" \
      > "${unit_root}/stdout.json" 2> "${unit_root}/stderr.txt" &
  pids+=("$!")
done
status=0
for pid in "${pids[@]}"; do
  if ! wait "${pid}"; then status=1; fi
done
if [[ "${status}" -ne 0 ]]; then
  printf 'At least one corrected convergence diagnostic unit failed.\n' >&2
  exit 1
fi

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

root = Path(sys.argv[1]).resolve(strict=True)
expected = {(1, 1701)}
observed = set()
for path in sorted((root / "units").glob("outer_*--seed_*")):
    outer = int(path.name.split("--")[0].removeprefix("outer_"))
    seed = int(path.name.split("--")[1].removeprefix("seed_"))
    receipt = json.loads((path / "fit/receipt.json").read_text(encoding="utf-8"))
    if (
        receipt.get("status") != "passed_timing_predictions_unscored"
        or receipt.get("outer_fold") != outer
        or receipt.get("seed") != seed
        or receipt.get("exact_production_grids_used") is not True
        or receipt.get("prediction_files") != 11
        or receipt.get("elastic_logistic_solver", {}).get("previous_max_iter") != 20000
        or receipt.get("elastic_logistic_solver", {}).get("intermediate_diagnostic_max_iter") != 100000
        or receipt.get("elastic_logistic_solver", {}).get("max_iter") != 500000
        or receipt.get("outer_test_outcomes_read") is not False
        or receipt.get("outer_test_metrics_calculated") is not False
    ):
        raise SystemExit(f"corrected convergence receipt differs: {path}")
    if len(list((path / "fit/predictions").glob("*.tsv"))) != 11:
        raise SystemExit(f"corrected convergence predictions differ: {path}")
    observed.add((outer, seed))
if observed != expected:
    raise SystemExit("corrected convergence unit roster differs")
PY
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/fit_gse267145_histology_baselines.py \
  tests/unit/test_fit_gse267145_histology_baselines.py \
  slurm/diagnose_gse267145_histology_convergence_retry_cpu.sbatch \
  "${MOLECULAR}/ARTIFACTS.json" \
  "${FOLDS}/ARTIFACTS.json" \
  "${FIT_VIEWS}/ARTIFACTS.json" \
  "${SURFACE}" \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    root,
    {
        "artifact_class": "gse267145_histology_convergence_retry_diagnostic",
        "logical_units": 1,
        "elastic_logistic_max_iter": 500000,
        "metrics_calculated": False,
        "outcomes_read": False,
        "production_authorized": False,
        "status": "passed_unscored_diagnostic",
    },
)
verify_frozen_tree(root)
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
sha256sum "${TARGET}/ARTIFACTS.json"
