#!/usr/bin/env bash
#SBATCH --job-name=masld-enformer-track-map
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
INSPECTION="${ROOT}/executions/enformer-crested-safe-keras-inspection-21065759"
SOURCES="${ROOT}/executions/borzoi-enformer-admission-sources-21064438"
TARGET="${ROOT}/executions/enformer-crested-track-crosswalk-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

cleanup() {
  if [[ -d "${STAGE}" ]]; then
    rm -rf -- "${STAGE}"
  fi
}
trap cleanup EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Enformer track crosswalk: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

module load python/3.11.5-GCCcore-13.2.0
python "${ROOT}/tests/unit/enformer_crested_test_track_crosswalk.py" \
  > "${STAGE}/validation/unittest.stdout.txt" \
  2> "${STAGE}/validation/unittest.stderr.txt"
python "${ROOT}/scripts/enformer_crested_build_track_crosswalk.py" \
  --labels "${INSPECTION}/inspection/output_classes.txt" \
  --targets "${SOURCES}/sources/enformer_targets_human.txt" \
  --output "${STAGE}/crosswalk" \
  > "${STAGE}/validation/crosswalk.stdout.txt"
python --version > "${STAGE}/validation/python.version.txt"
pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${INSPECTION}/ARTIFACTS.json" \
  "${SOURCES}/ARTIFACTS.json" \
  "${ROOT}/scripts/enformer_crested_build_track_crosswalk.py" \
  "${ROOT}/tests/unit/enformer_crested_test_track_crosswalk.py" \
  "${ROOT}/slurm/enformer_crested_build_track_crosswalk.sbatch" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "crosswalk/receipt.json").read_text())
freeze_tree(
    root,
    {
        "artifact_class": "restricted_sequence_model_output_crosswalk",
        "model_id": receipt["model_id"],
        "canonical_target_rows": receipt["canonical_target_rows"],
        "position_identity": receipt["position_by_position_description_identity"],
        "native_parity_established": False,
        "outcomes_read": False,
        "open_champion_eligible": False,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Enformer CREsted track crosswalk: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
