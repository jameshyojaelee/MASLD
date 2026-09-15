#!/usr/bin/env bash
#SBATCH --job-name=masld-corgi-smoke-tiles
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TILES="${ROOT}/executions/corgi-outcome-aligned-tiles-21068311"
TARGET="${ROOT}/executions/corgi-tile-smoke-subset-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Corgi smoke subset artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
printf '%s  %s\n' \
  1cdf94aa4c37ead139125ffc3d8fee358605a2394a206b1f3e79190ecd25645c \
  "${TILES}/ARTIFACTS.json" | sha256sum --check --strict
install -d -m 0750 "${STAGE}/validation"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}" python -m unittest tests.unit.test_build_corgi_tile_smoke_subset \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}" python "${ROOT}/scripts/build_corgi_tile_smoke_subset.py" \
  --tile-contract "${TILES}" --tiles-per-fold 128 --output "${STAGE}/subset" \
  > "${STAGE}/validation/build.stdout.json" \
  2> "${STAGE}/validation/build.stderr.txt"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/build_corgi_tile_smoke_subset.py" \
  "${ROOT}/tests/unit/test_build_corgi_tile_smoke_subset.py" \
  "${ROOT}/slurm/build_corgi_tile_smoke_subset.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads((root / "subset" / "summary.json").read_text())
if (
    summary["status"] != "pass_outcome_free_smoke_subset"
    or summary["tiles"] != 640
    or summary["outcomes_or_labels_read"]
    or summary["model_predictions_read"]
):
    raise SystemExit("Corgi smoke subset receipt differs")
freeze_tree(root, {
    "artifact_class": "corgi_outcome_aligned_tile_smoke_subset",
    "dataset_id": "gse296875",
    "tiles": summary["tiles"],
    "windows": summary["windows"],
    "outcomes_read": False,
    "predictions_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Corgi tile smoke subset: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
