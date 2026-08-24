#!/usr/bin/env bash
#SBATCH --job-name=masld-borzoi-enformer-test
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/borzoi-enformer-admission-test-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Borzoi/Enformer admission test artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="${ROOT}" python -m unittest -v \
  tests.unit.borzoi_enformer_test_fixture \
  > "${STAGE}/unittest.stdout.txt" \
  2> "${STAGE}/unittest.stderr.txt"
python --version > "${STAGE}/python.version.txt" 2>&1
python -m pip freeze > "${STAGE}/pip-freeze.txt"
sha256sum \
  "${ROOT}/config/borzoi_enformer_admission.json" \
  "${ROOT}/scripts/borzoi_enformer_build_fixture.py" \
  "${ROOT}/scripts/borzoi_enformer_validate_sources.py" \
  "${ROOT}/tests/unit/borzoi_enformer_test_fixture.py" \
  "${ROOT}/slurm/borzoi_enformer_test_admission.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

freeze_tree(
    Path(sys.argv[1]).resolve(strict=True),
    {
        "artifact_class": "borzoi_enformer_admission_test",
        "models": ["borzoi_ensemble", "enformer"],
        "checkpoint_bytes_loaded": False,
        "observed_outcomes_loaded": False,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'Borzoi/Enformer admission test artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
