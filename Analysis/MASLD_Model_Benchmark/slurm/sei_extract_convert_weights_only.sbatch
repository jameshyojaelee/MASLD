#!/usr/bin/env bash
#SBATCH --job-name=masld-sei-weights-only
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python
SOURCE="${ROOT}/executions/alphagenome-sei-admission-sources-21064512"
SOURCE_SHA=84ae3470dca123c20904e315ef9926d5a139fb9fe6d99c51f6aaf694ef6e6785
ARCHIVE="${ROOT}/executions/sei-safe-native-archive-21064521"
ARCHIVE_SHA=beff13e916cf12eb6db2e29cb301b3151cc519048293a2ada2a66f6fa3691d2d
RUNTIME="${ROOT}/executions/environments/dna-lm-transformers-4.38.1-21064683"
RUNTIME_SHA=876c072889324661df0500e2f1eb8d09a8f1e3054deddb153d8296b069845a63
CHECKPOINT_SHA=9d4771d5363e4955eec2a36567fedac986191f64c7874c85bfb9e203c7f60783
CHECKPOINT_SIZE=3559928168
TARGET="${ROOT}/executions/sei-weights-only-conversion-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed Sei weights-only stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Sei weights-only artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/validation" \
  "${STAGE}/offline-guard" \
  "${STAGE}/hf-home"
for binding in \
  "${SOURCE_SHA} ${SOURCE}" \
  "${ARCHIVE_SHA} ${ARCHIVE}" \
  "${RUNTIME_SHA} ${RUNTIME}"; do
  read -r expected directory <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${directory}/ARTIFACTS.json" \
    | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${SOURCE}" "${ARCHIVE}" "${RUNTIME}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
summary = json.loads(
    (Path(sys.argv[2]) / "inventory" / "admission_summary.json").read_text()
)
if (
    summary["status"] != "pass"
    or summary["checkpoint_deserialized"]
    or summary["next_allowed_action"]
    != "isolated_no_network_tensor_only_checkpoint_inspection"
):
    raise SystemExit("Sei archive does not permit tensor-only inspection")
print("Sei frozen inputs and next-action gate verified")
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
"${BASE}" "${ROOT}/scripts/sei_extract_native_members.py" \
  --archive "${ARCHIVE}/checkpoint/sei_model.tar.gz" \
  --admission-summary "${ARCHIVE}/inventory/admission_summary.json" \
  --output "${STAGE}/native" \
  > "${STAGE}/validation/extraction.stdout.txt" \
  2> "${STAGE}/validation/extraction.stderr.txt"
"${BASE}" "${ROOT}/scripts/torch_legacy_checkpoint_safe_inventory.py" \
  --input "${STAGE}/native/sei.pth" \
  --output "${STAGE}/validation/checkpoint_container_inventory.json" \
  --expected-sha256 "${CHECKPOINT_SHA}" \
  --expected-size "${CHECKPOINT_SIZE}" \
  > "${STAGE}/validation/container.stdout.txt" \
  2> "${STAGE}/validation/container.stderr.txt"
"${BASE}" "${ROOT}/scripts/dna_lm_convert_weights_only.py" \
  --source "${STAGE}/native/sei.pth" \
  --output "${STAGE}/converted" \
  --expected-sha256 "${CHECKPOINT_SHA}" \
  --expected-size "${CHECKPOINT_SIZE}" \
  --model-id sei \
  --source-container PyTorch_legacy_pickle_stream \
  --disable-mmap \
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
  --input "${STAGE}/converted/sei.model.safetensors" \
  --output "${STAGE}/converted/safetensors_inventory.json" \
  --expected-sha256 "${CONVERTED_SHA}" \
  --expected-size "${CONVERTED_SIZE}" \
  > "${STAGE}/validation/inventory.stdout.txt" \
  2> "${STAGE}/validation/inventory.stderr.txt"
rm "${STAGE}/native/sei.pth"
"${BASE}" -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/sei_extract_native_members.py" \
  "${ROOT}/scripts/torch_legacy_checkpoint_safe_inventory.py" \
  "${ROOT}/scripts/dna_lm_convert_weights_only.py" \
  "${ROOT}/scripts/safetensors_safe_inventory.py" \
  "${ROOT}/slurm/sei_extract_convert_weights_only.sbatch" \
  "${SOURCE}/ARTIFACTS.json" \
  "${ARCHIVE}/ARTIFACTS.json" \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
extraction = json.loads((root / "native" / "extraction_receipt.json").read_text())
conversion = json.loads((root / "converted" / "conversion_receipt.json").read_text())
inventory = json.loads((root / "converted" / "safetensors_inventory.json").read_text())
guard = json.loads((root / "validation" / "offline_guard.json").read_text())
if (
    extraction["status"] != "pass"
    or extraction["projection_shape"] != [61, 21907]
    or extraction["projection_basis_count"] != 61
    or extraction["sequence_class_count"] != 40
    or not extraction["projection_to_sequence_class_mapping_unresolved"]
    or extraction["checkpoint_deserialized"]
    or extraction["model_forward_executed"]
    or conversion["status"] != "pass"
    or conversion["model_id"] != "sei"
    or conversion["source_container"] != "PyTorch_legacy_pickle_stream"
    or conversion["loader"] != "torch.load_weights_only_true_mmap_false"
    or conversion["explicit_weights_only_safe_globals"]
    or conversion["custom_model_code_imported"]
    or conversion["model_instantiated"]
    or conversion["model_forward_executed"]
    or conversion["observed_outcomes_loaded"]
    or conversion["sealed_outcomes_loaded"]
    or inventory["tensor_data_loaded"]
    or not guard["python_INET_socket_blocked"]
    or (root / "native" / "sei.pth").exists()
):
    raise SystemExit("Sei weights-only conversion receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "sei_weights_only_conversion",
        "model_id": "sei",
        "source_artifacts_sha256": "84ae3470dca123c20904e315ef9926d5a139fb9fe6d99c51f6aaf694ef6e6785",
        "archive_artifacts_sha256": "beff13e916cf12eb6db2e29cb301b3151cc519048293a2ada2a66f6fa3691d2d",
        "runtime_artifacts_sha256": "876c072889324661df0500e2f1eb8d09a8f1e3054deddb153d8296b069845a63",
        "loader": "torch_load_weights_only_true_mmap_false_legacy_stream",
        "terms": "official_code_research_only_restricted_comparator",
        "outcome_firewall": "no_observed_or_sealed_outcomes",
        "custom_model_code_imported": False,
        "model_forward_executed": False,
        "open_champion_eligible": False,
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Sei weights-only conversion: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
