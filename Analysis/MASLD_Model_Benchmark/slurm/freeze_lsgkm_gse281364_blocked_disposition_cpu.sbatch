#!/usr/bin/env bash
#SBATCH --job-name=model-check-284
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONFIG="${ROOT}/config/lsgkm_gse281364_blocked_disposition.json"
SCRIPT="${ROOT}/scripts/freeze_lsgkm_gse281364_blocked_disposition.py"
TARGET="${ROOT}/executions/model-check-284-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite blocked disposition: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

printf '%s  %s\n' 15a7aaadca2064c0100ae6d07423f32692f91cf9b7ee5793722b94ac9915d375 \
  "${CONFIG}" | sha256sum --check --strict
printf '%s  %s\n' 560a213e05662705d70d61147f667d6a719a99625871b284bb6497a1e559c2a5 \
  "${SCRIPT}" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
python -c \
  'from pathlib import Path; import sys; compile(Path(sys.argv[1]).read_text(encoding="utf-8"), sys.argv[1], "exec")' \
  "${SCRIPT}"
python "${SCRIPT}" \
  --project-root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/terminal_disposition.json" \
  > "${STAGE}/audit.stdout.json" \
  2> "${STAGE}/audit.stderr.txt"
python -m pip freeze --all > "${STAGE}/pip-freeze.txt"
sha256sum "${SCRIPT}" "${CONFIG}" > "${STAGE}/source.sha256"
module purge

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
disposition = json.loads((root / "terminal_disposition.json").read_text(encoding="utf-8"))
failure = disposition["first_failed_fit"]
if (
    disposition["status"] != "terminally_blocked_300bp_matched_negative_design"
    or failure["gate"]["matched_pairs"] != 2097
    or failure["gate"]["positive_coverage"] != 0.028185862713208512
    or disposition["input_fastas_written"] != 0
    or disposition["production_fits_executed"] != 0
    or disposition["production_predictions_generated"] != 0
    or disposition["outcomes_read"]
    or disposition["reporter_counts_read"]
    or disposition["sealed_assets_read"]
):
    raise SystemExit("terminal LS-GKM disposition differs")
freeze_tree(root, {
    "artifact_class": "lsgkm_gse281364_terminal_blocked_disposition",
    "planned_shared_fit_count": 25,
    "evaluated_fit_units": 1,
    "matched_pairs_first_failed_fit": 2097,
    "input_fastas_written": 0,
    "production_fits_executed": 0,
    "production_predictions_generated": 0,
    "outcomes_read": False,
    "reporter_counts_read": False,
    "sealed_assets_read": False,
    "status": "terminally_blocked_300bp_matched_negative_design",
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Terminal blocked disposition: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
