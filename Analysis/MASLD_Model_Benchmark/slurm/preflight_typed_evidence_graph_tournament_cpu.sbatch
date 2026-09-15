#!/usr/bin/env bash
#SBATCH --job-name=model-check-217
#SBATCH --partition=cpu
#SBATCH --account=nslab
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
CONFIG="${ROOT}/config/typed_evidence_graph_tournament_preflight.json"
TARGET="${ROOT}/executions/model-check-217-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite graph preflight: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
export PYTHONPATH="${ROOT}:${ROOT}/src"

"${PYTHON}" -m unittest -v \
  tests.unit.test_preflight_typed_evidence_graph_tournament \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

"${PYTHON}" "${ROOT}/scripts/preflight_typed_evidence_graph_tournament.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/preflight" \
  > "${STAGE}/validation/preflight.stdout.json" \
  2> "${STAGE}/validation/preflight.stderr.txt"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/preflight_typed_evidence_graph_tournament.py" \
  "${ROOT}/tests/unit/test_preflight_typed_evidence_graph_tournament.py" \
  "${ROOT}/slurm/preflight_typed_evidence_graph_tournament_cpu.sbatch" \
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
if (
    receipt["status"]
    != "blocked_no_executable_wholly_removed_source_task"
    or receipt["strongest_currently_executable_task"] is not None
    or receipt["model_census_rows"] != 12
    or receipt["mandatory_baseline_rows"] != 5
    or receipt["candidate_source_rows"] != 3
    or receipt["activation_gate_rows"] != 17
    or receipt["required_gates_unresolved"] != 16
    or receipt["sealed_graph_source_registered"]
    or receipt["graph_snapshot_bound"]
    or receipt["wholly_removed_source_bound"]
    or receipt["explicit_tested_negative_denominator_bound"]
    or receipt["training_or_scoring_allowed"]
    or receipt["champion_claim_allowed"]
    or receipt["raw_or_summary_labels_read"]
    or receipt["sealed_assets_read"]
    or receipt["model_import_fit_predict_or_score"]
    or receipt["network_download_or_install"]
    or receipt["bundled_cpu_validation_ask"]["submitted"]
):
    raise SystemExit("typed-evidence graph preflight receipt differs")

freeze_tree(
    stage,
    {
        "artifact_class": "typed_evidence_graph_tournament_preflight",
        "model_count": 12,
        "mandatory_baseline_count": 5,
        "candidate_source_count": 3,
        "activation_gate_count": 17,
        "training_or_scoring_allowed": False,
        "sealed_assets_read": False,
        "status": "blocked_no_executable_wholly_removed_source_task",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Graph preflight: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
