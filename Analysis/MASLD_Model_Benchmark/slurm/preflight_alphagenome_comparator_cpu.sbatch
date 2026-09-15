#!/usr/bin/env bash
#SBATCH --job-name=model-check-215
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
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
CONFIG="${ROOT}/config/alphagenome_comparator_preflight.json"
TARGET="${ROOT}/executions/model-check-215-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite comparator preflight: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
export PYTHONPATH="${ROOT}:${ROOT}/src"

"${PYTHON}" -m unittest -v \
  tests.unit.test_preflight_alphagenome_comparator \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

"${PYTHON}" "${ROOT}/scripts/preflight_alphagenome_comparator.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/preflight" \
  > "${STAGE}/validation/preflight.stdout.json" \
  2> "${STAGE}/validation/preflight.stderr.txt"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/preflight_alphagenome_comparator.py" \
  "${ROOT}/tests/unit/test_preflight_alphagenome_comparator.py" \
  "${ROOT}/slurm/preflight_alphagenome_comparator_cpu.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree


stage = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads(
    (stage / "preflight/receipt.json").read_text(encoding="utf-8")
)
metadata = receipt["metadata_audit"]
local = receipt["local_checkpoint_route"]
api = receipt["api_route"]
if (
    receipt["status"]
    != "terminal_current_release_restricted_and_execution_blocked"
    or receipt["executable"]
    or not receipt["terminal_current_release"]
    or metadata["encoded_output_slots"] != 6126
    or metadata["padding_slots"] != 563
    or metadata["nonpadding_metadata_rows"] != 5563
    or metadata["splice_junction_metadata_rows"] != 367
    or metadata["active_tracks_with_second_splice_junction_strands"] != 5930
    or receipt["relevant_track_roster_rows"] != 93
    or len(receipt["blockers"]) != 21
    or local["checkpoint_downloaded"]
    or local["checkpoint_deserialized"]
    or local["model_forward_executed"]
    or api["api_key_inspected"]
    or api["api_connection_attempted"]
    or api["metadata_endpoint_called"]
    or api["prediction_endpoint_called"]
    or api["immutable_server_build_identity_available"]
    or receipt["outcomes_read"]
    or receipt["sealed_assets_read"]
    or receipt["metrics_calculated"]
    or receipt["open_champion_eligible"]
    or receipt["conditional_open_model_training_allowed"]
):
    raise SystemExit("restricted comparator receipt differs")

freeze_tree(
    stage,
    {
        "artifact_class": "restricted_sequence_comparator_preflight",
        "model_id": "alphagenome",
        "executable": False,
        "model_forward_executed": False,
        "api_connection_attempted": False,
        "outcomes_read": False,
        "sealed_assets_read": False,
        "open_champion_eligible": False,
        "status": "terminal_current_release_restricted_and_execution_blocked",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Comparator preflight: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
