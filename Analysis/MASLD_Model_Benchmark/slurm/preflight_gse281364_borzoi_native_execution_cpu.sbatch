#!/usr/bin/env bash
#SBATCH --job-name=model-check-214
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
CONFIG="${ROOT}/config/gse281364_borzoi_native_execution.json"
TARGET="${ROOT}/executions/borzoi-native-execution-preflight-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite native execution preflight: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m unittest -v \
  tests.unit.test_preflight_gse281364_borzoi_native_execution \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${PYTHON}" "${ROOT}/scripts/preflight_gse281364_borzoi_native_execution.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/receipt.json" \
  --mode audit \
  > "${STAGE}/validation/preflight.stdout.json" \
  2> "${STAGE}/validation/preflight.stderr.txt"
"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/preflight_gse281364_borzoi_native_execution.py" \
  "${ROOT}/scripts/run_gse281364_borzoi_native_predictions.py" \
  "${ROOT}/tests/unit/test_preflight_gse281364_borzoi_native_execution.py" \
  "${ROOT}/slurm/preflight_gse281364_borzoi_native_execution_cpu.sbatch" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
import json
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((stage / "receipt.json").read_text(encoding="utf-8"))
gates = receipt["execution_gates"]
expected = {
    "weight_terms",
    "official_checkpoints",
    "native_runtime",
    "reference_parity",
    "native_parity",
}
if (
    receipt["status"]
    != "blocked_missing_terms_weights_runtime_reference_parity_native_parity"
    or receipt["executable"]
    or set(gates) != expected
    or any(value["passed"] for value in gates.values())
    or len(receipt["blockers"]) != 11
    or receipt["local_port_observation"]["eligible_as_native_runtime"]
    or receipt["checkpoint_bytes_opened"]
    or receipt["model_modules_imported"]
    or receipt["model_forward_executed"]
    or receipt["reporter_counts_read"]
    or receipt["reporter_outcomes_read"]
    or receipt["sealed_assets_read"]
    or receipt["metrics_calculated"]
):
    raise SystemExit("native Borzoi blocked preflight receipt differs")
freeze_tree(
    stage,
    {
        "artifact_class": "gse281364_borzoi_native_execution_preflight",
        "dataset_id": "gse281364",
        "model_id": "borzoi_ensemble",
        "executable": False,
        "checkpoint_bytes_opened": False,
        "model_forward_executed": False,
        "outcomes_read": False,
        "status": "blocked_missing_terms_weights_runtime_reference_parity_native_parity",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Native execution preflight: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
