#!/usr/bin/env bash
#SBATCH --job-name=model-check-077
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
TARGET=${ROOT}/executions/model-check-077-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failure() {
  exit_code=$?
  if [[ "${exit_code}" -ne 0 && -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
  exit "${exit_code}"
}
trap preserve_failure EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE267145 production-bundle validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/validation/input.sha256"
sha256sum --check --strict "${STAGE}/validation/input.sha256"
bash -n slurm/run_gse267145_histology_production_bundle_cpu.sbatch
PYTHONPATH="${ROOT}:${ROOT}/src" "${RUNTIME}/env/bin/python" -m unittest -v \
  tests.unit.test_fit_gse267145_histology_baselines \
  tests.unit.test_aggregate_gse267145_histology_production_predictions \
  tests.unit.test_gse267145_histology_production_bundle_contract \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${RUNTIME}/env/bin/python" -m pip freeze --all \
  > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/fit_gse267145_histology_baselines.py \
  scripts/run_gse267145_histology_production_bundle.py \
  scripts/aggregate_gse267145_histology_production_predictions.py \
  tests/unit/test_fit_gse267145_histology_baselines.py \
  tests/unit/test_aggregate_gse267145_histology_production_predictions.py \
  tests/unit/test_gse267145_histology_production_bundle_contract.py \
  slurm/run_gse267145_histology_production_bundle_cpu.sbatch \
  slurm/test_gse267145_histology_production_bundle_cpu.sbatch \
  > "${STAGE}/validation/source.sha256"
PYTHONPATH="${ROOT}:${ROOT}/src" "${RUNTIME}/env/bin/python" - "${STAGE}" <<'PY'
from pathlib import Path
import re
import sys

from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
tests = (stage / "validation/tests.stderr.txt").read_text(encoding="utf-8")
match = re.search(r"Ran (\d+) tests", tests)
if match is None or int(match.group(1)) != 19 or "OK" not in tests:
    raise SystemExit("GSE267145 production-bundle tests did not pass")
freeze_tree(
    stage,
    {
        "artifact_class": "gse267145_histology_production_bundle_validation",
        "tests_passed": int(match.group(1)),
        "production_submitted_by_validation": False,
        "metrics_calculated": False,
        "outcomes_read": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
