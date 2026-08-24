#!/usr/bin/env bash
#SBATCH --job-name=masld-enformer-port-inspect
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
ADMISSION="${ROOT}/executions/enformer-crested-restricted-port-admission-21065579"
TARGET="${ROOT}/executions/enformer-crested-safe-keras-inspection-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

cleanup() {
  if [[ -d "${STAGE}" ]]; then
    rm -rf -- "${STAGE}"
  fi
}
trap cleanup EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Enformer Keras inspection: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/scratch" "${STAGE}/validation"

module load python/3.11.5-GCCcore-13.2.0
python "${ROOT}/tests/unit/enformer_crested_test_inspect_keras.py" \
  > "${STAGE}/validation/unittest.stdout.txt" \
  2> "${STAGE}/validation/unittest.stderr.txt"
python "${ROOT}/scripts/enformer_crested_inspect_keras.py" \
  --archive "${ADMISSION}/downloads/enformer_human.tar.gz" \
  --output "${STAGE}/inspection" \
  --scratch "${STAGE}/scratch" \
  > "${STAGE}/validation/inspection.stdout.txt"
rmdir "${STAGE}/scratch"
python --version > "${STAGE}/validation/python.version.txt"
pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${ADMISSION}/ARTIFACTS.json" \
  "${ROOT}/scripts/enformer_crested_inspect_keras.py" \
  "${ROOT}/tests/unit/enformer_crested_test_inspect_keras.py" \
  "${ROOT}/slurm/enformer_crested_inspect_keras.sbatch" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "inspection/receipt.json").read_text())
freeze_tree(
    root,
    {
        "artifact_class": "restricted_keras_checkpoint_safe_inspection",
        "model_id": receipt["model_id"],
        "outer_archive_sha256": receipt["outer_archive_sha256"],
        "keras_deserialization_performed": False,
        "model_forward_executed": False,
        "outcomes_read": False,
        "open_champion_eligible": False,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Enformer CREsted safe Keras inspection: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
