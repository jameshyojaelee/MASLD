#!/usr/bin/env bash
#SBATCH --job-name=masld-enformer-eval-preflight
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON="${PROJECT}/.agent/scanpy_env/bin/python"
TARGET="${ROOT}/executions/enformer-evaluator-preflight-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed Enformer evaluator preflight: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Enformer evaluator preflight: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

PYTHONPATH="${ROOT}" "${PYTHON}" -m unittest -v \
  tests.unit.test_evaluate_gse281364_enformer \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}" "${PYTHON}" \
  "${ROOT}/scripts/evaluate_gse281364_enformer.py" --help \
  > "${STAGE}/validation/cli-help.stdout.txt" \
  2> "${STAGE}/validation/cli-help.stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/evaluate_gse281364_enformer.py" \
  "${ROOT}/tests/unit/test_evaluate_gse281364_enformer.py" \
  "${ROOT}/slurm/preflight_gse281364_enformer_evaluator_cpu.sbatch" \
  "${ROOT}/slurm/evaluate_gse281364_enformer_cpu.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
tests = (root / "validation/tests.stderr.txt").read_text(encoding="utf-8")
help_text = (root / "validation/cli-help.stdout.txt").read_text(encoding="utf-8")
if (
    "Ran 7 tests" not in tests
    or "OK" not in tests
    or "--predictions" not in help_text
    or "--outcomes" not in help_text
    or (root / "validation/cli-help.stderr.txt").stat().st_size != 0
):
    raise SystemExit("Enformer evaluator source/test preflight differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse281364_enformer_evaluator_outcome_free_preflight",
        "dataset_id": "gse281364",
        "model_id": "enformer_crested_restricted_port",
        "synthetic_tests": 7,
        "cli_help_passed": True,
        "project_outcomes_read": False,
        "prediction_artifact_read": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Enformer evaluator outcome-free preflight: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
