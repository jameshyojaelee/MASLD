#!/usr/bin/env bash
#SBATCH --job-name=masld-evo-weights-only
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
INPUT="${ROOT}/executions/dna-lm-safe-pickle-checkpoints-21064576"
INPUT_SHA=0ab85bab69cdc5e27aa2321c2664050f84f30a4fa077fc89a380673eecec6edb
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
BASE_SHA=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
TORCH_SHA=d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3
CONVERSION="${ROOT}/executions/environments/dna-lm-torch2.6-conversion-21064945"
CONVERSION_SHA=5fed22c975200495ea53b38b09f5286cf13a9230a138cb72a0c758c85d5cb9d4
CHECKPOINT="${INPUT}/checkpoints/evo.evo2_7b.pt"
CHECKPOINT_SHA=c66645929dc1b9c631f5be656da8726f38946315dc9167000a615dd626fcecf4
CHECKPOINT_SIZE=13766621200
TARGET="${ROOT}/executions/dna-lm-weights-only-evo-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed Evo weights-only stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Evo weights-only artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/validation" \
  "${STAGE}/offline-guard" \
  "${STAGE}/hf-home"
for binding in \
  "${INPUT_SHA} ${INPUT}" \
  "${BASE_SHA} ${BASE}" \
  "${TORCH_SHA} ${TORCH}" \
  "${CONVERSION_SHA} ${CONVERSION}"; do
  read -r expected directory <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${directory}/ARTIFACTS.json" \
    | sha256sum --check --strict
done
printf '%s  %s\n' "${CHECKPOINT_SHA}" "${CHECKPOINT}" \
  | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${INPUT}" "${BASE}" "${TORCH}" "${CONVERSION}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("Evo weights-only frozen inputs verified")
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
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${STAGE}/offline-guard:${CONVERSION}/site-packages:${TORCH}/site-packages:${BASE}/site-packages"
python - <<'PY' > "${STAGE}/validation/offline_guard.json"
import json
import socket
import torch

blocked = False
try:
    socket.socket(socket.AF_INET, socket.SOCK_STREAM)
except RuntimeError:
    blocked = True
receipt = {
    "python_INET_socket_blocked": blocked,
    "torch": str(torch.__version__),
    "torch_cuda": torch.version.cuda,
    "weights_only_safe_globals_api": hasattr(
        torch.serialization, "add_safe_globals"
    ),
}
print(json.dumps(receipt, sort_keys=True))
if (
    not blocked
    or receipt["torch"] != "2.6.0+cu124"
    or receipt["torch_cuda"] != "12.4"
    or not receipt["weights_only_safe_globals_api"]
):
    raise SystemExit(f"Evo conversion runtime or network guard differs: {receipt!r}")
PY
python "${ROOT}/scripts/dna_lm_convert_weights_only.py" \
  --source "${CHECKPOINT}" \
  --output "${STAGE}/converted" \
  --expected-sha256 "${CHECKPOINT_SHA}" \
  --expected-size "${CHECKPOINT_SIZE}" \
  --model-id evo \
  > "${STAGE}/validation/conversion.stdout.txt" \
  2> "${STAGE}/validation/conversion.stderr.txt"
read -r CONVERTED_SHA CONVERTED_SIZE < <(
  python - "${STAGE}/converted/conversion_receipt.json" <<'PY'
import json
import sys

receipt = json.load(open(sys.argv[1]))
print(receipt["converted_safetensors_sha256"], receipt["converted_safetensors_size_bytes"])
PY
)
python "${ROOT}/scripts/safetensors_safe_inventory.py" \
  --input "${STAGE}/converted/evo.model.safetensors" \
  --output "${STAGE}/converted/safetensors_inventory.json" \
  --expected-sha256 "${CONVERTED_SHA}" \
  --expected-size "${CONVERTED_SIZE}" \
  > "${STAGE}/validation/inventory.stdout.txt" \
  2> "${STAGE}/validation/inventory.stderr.txt"
python -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/dna_lm_convert_weights_only.py" \
  "${ROOT}/scripts/safetensors_safe_inventory.py" \
  "${ROOT}/slurm/dna_lm_convert_weights_only_evo_torch26.sbatch" \
  "${INPUT}/ARTIFACTS.json" \
  "${BASE}/ARTIFACTS.json" \
  "${TORCH}/ARTIFACTS.json" \
  "${CONVERSION}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
module purge

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "converted" / "conversion_receipt.json").read_text())
inventory = json.loads((root / "converted" / "safetensors_inventory.json").read_text())
guard = json.loads((root / "validation" / "offline_guard.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["model_id"] != "evo"
    or receipt["explicit_weights_only_safe_globals"] != ["_io.BytesIO"]
    or receipt["custom_model_code_imported"]
    or receipt["model_instantiated"]
    or receipt["model_forward_executed"]
    or receipt["observed_outcomes_loaded"]
    or receipt["sealed_outcomes_loaded"]
    or inventory["tensor_data_loaded"]
    or not guard["python_INET_socket_blocked"]
    or not guard["weights_only_safe_globals_api"]
):
    raise SystemExit("Evo weights-only conversion receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "dna_lm_weights_only_conversion",
        "model_id": "evo",
        "input_artifacts_sha256": "0ab85bab69cdc5e27aa2321c2664050f84f30a4fa077fc89a380673eecec6edb",
        "base_artifacts_sha256": "30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af",
        "torch_artifacts_sha256": "d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3",
        "conversion_artifacts_sha256": "5fed22c975200495ea53b38b09f5286cf13a9230a138cb72a0c758c85d5cb9d4",
        "loader": "torch2.6_load_weights_only_true_mmap_true",
        "explicit_safe_globals": ["_io.BytesIO"],
        "outcome_firewall": "no_observed_or_sealed_outcomes",
        "custom_model_code_imported": False,
        "model_forward_executed": False,
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Evo weights-only conversion: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
