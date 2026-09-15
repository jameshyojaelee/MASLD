#!/usr/bin/env bash
#SBATCH --job-name=corgi-regional-eval-test
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/corgi-regional-evaluator-test-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite evaluator test artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' \
  8008d5a8b8803b61890db4823d4b27be4c472f44881ca3788ce728f01be8f731 \
  "${ROOT}/scripts/evaluate_corgi_regional_smoke.py" | sha256sum --check --strict
printf '%s  %s\n' \
  f2725c66fae08cb4e2002008441a36e30b8591abc6cdc11731fff1fa3893fa16 \
  "${ROOT}/tests/unit/test_evaluate_corgi_regional_smoke.py" | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
cd "${ROOT}"
PYTHONPATH="src:." python -m unittest -v tests.unit.test_evaluate_corgi_regional_smoke \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  scripts/evaluate_corgi_regional_smoke.py \
  tests/unit/test_evaluate_corgi_regional_smoke.py \
  slurm/test_evaluate_corgi_regional_smoke.sbatch \
  > "${STAGE}/source.sha256"
micromamba deactivate
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    root,
    {
        "artifact_class": "corgi_regional_smoke_evaluator_unit_test",
        "tests": 5,
        "outcomes_read": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Corgi regional evaluator tests: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
