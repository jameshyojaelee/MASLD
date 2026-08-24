#!/usr/bin/env bash
set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

if [[ "$#" -ne 5 ]]; then
  printf 'Usage: %s MODEL_ID FILENAME SHA256 SIZE WRAPPER\n' "$0" >&2
  exit 2
fi
MODEL_ID="$1"
FILENAME="$2"
CHECKPOINT_SHA="$3"
CHECKPOINT_SIZE="$4"
WRAPPER="$5"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python
INPUT="${ROOT}/executions/dna-lm-safe-pickle-checkpoints-21064576"
INPUT_SHA=0ab85bab69cdc5e27aa2321c2664050f84f30a4fa077fc89a380673eecec6edb
RUNTIME="${ROOT}/executions/environments/dna-lm-transformers-4.38.1-21064683"
RUNTIME_SHA=876c072889324661df0500e2f1eb8d09a8f1e3054deddb153d8296b069845a63
TARGET="${ROOT}/executions/dna-lm-weights-only-${MODEL_ID}-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed weights-only stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite weights-only artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/validation" \
  "${STAGE}/offline-guard" \
  "${STAGE}/hf-home"
printf '%s  %s\n' "${INPUT_SHA}" "${INPUT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${RUNTIME_SHA}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${CHECKPOINT_SHA}" "${INPUT}/checkpoints/${FILENAME}" \
  | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${INPUT}" "${RUNTIME}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("weights-only frozen inputs verified")
PY
module purge

cat > "${STAGE}/offline-guard/sitecustomize.py" <<'PY'
import socket

_socket = socket.socket

class GuardedSocket(_socket):
    def __new__(cls, family=-1, *args, **kwargs):
        if family in {socket.AF_INET, socket.AF_INET6}:
            raise RuntimeError("network disabled for checkpoint conversion")
        return super().__new__(cls, family, *args, **kwargs)

socket.socket = GuardedSocket

def blocked_connection(*args, **kwargs):
    raise RuntimeError("network disabled for checkpoint conversion")

socket.create_connection = blocked_connection
PY
export HF_HOME="${STAGE}/hf-home"
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export PYTHONPATH="${STAGE}/offline-guard:${RUNTIME}/site-packages"
"${BASE}" - <<'PY' > "${STAGE}/validation/offline_guard.json"
import json
import socket

blocked = False
try:
    socket.socket(socket.AF_INET, socket.SOCK_STREAM)
except RuntimeError:
    blocked = True
if not blocked:
    raise SystemExit("network guard differs")
print(json.dumps({"python_INET_socket_blocked": blocked}, sort_keys=True))
PY

"${BASE}" "${ROOT}/scripts/dna_lm_convert_weights_only.py" \
  --source "${INPUT}/checkpoints/${FILENAME}" \
  --output "${STAGE}/converted" \
  --expected-sha256 "${CHECKPOINT_SHA}" \
  --expected-size "${CHECKPOINT_SIZE}" \
  --model-id "${MODEL_ID}" \
  > "${STAGE}/validation/conversion.stdout.txt" \
  2> "${STAGE}/validation/conversion.stderr.txt"
read -r CONVERTED_SHA CONVERTED_SIZE < <(
  "${BASE}" - "${STAGE}/converted/conversion_receipt.json" <<'PY'
import json
import sys

receipt = json.load(open(sys.argv[1]))
print(receipt["converted_safetensors_sha256"], receipt["converted_safetensors_size_bytes"])
PY
)
"${BASE}" "${ROOT}/scripts/safetensors_safe_inventory.py" \
  --input "${STAGE}/converted/${MODEL_ID}.model.safetensors" \
  --output "${STAGE}/converted/safetensors_inventory.json" \
  --expected-sha256 "${CONVERTED_SHA}" \
  --expected-size "${CONVERTED_SIZE}" \
  > "${STAGE}/validation/inventory.stdout.txt" \
  2> "${STAGE}/validation/inventory.stderr.txt"
"${BASE}" -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/dna_lm_convert_weights_only.py" \
  "${ROOT}/scripts/dna_lm_convert_weights_only_driver.sh" \
  "${ROOT}/scripts/safetensors_safe_inventory.py" \
  "${WRAPPER}" \
  "${INPUT}/ARTIFACTS.json" \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${MODEL_ID}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
model_id = sys.argv[2]
receipt = json.loads((root / "converted" / "conversion_receipt.json").read_text())
inventory = json.loads((root / "converted" / "safetensors_inventory.json").read_text())
guard = json.loads((root / "validation" / "offline_guard.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["model_id"] != model_id
    or receipt["custom_model_code_imported"]
    or receipt["model_instantiated"]
    or receipt["model_forward_executed"]
    or inventory["tensor_data_loaded"]
    or not guard["python_INET_socket_blocked"]
):
    raise SystemExit("weights-only conversion receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "dna_lm_weights_only_conversion",
        "model_id": model_id,
        "input_artifacts_sha256": "0ab85bab69cdc5e27aa2321c2664050f84f30a4fa077fc89a380673eecec6edb",
        "runtime_artifacts_sha256": "876c072889324661df0500e2f1eb8d09a8f1e3054deddb153d8296b069845a63",
        "loader": "torch.load_weights_only_true_mmap_true",
        "outcome_firewall": "no_observed_or_sealed_outcomes",
        "custom_model_code_imported": False,
        "model_forward_executed": False,
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'DNA-LM weights-only conversion: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
