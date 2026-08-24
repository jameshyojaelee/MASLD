#!/usr/bin/env bash
#SBATCH --job-name=masld-dna-lm-safe
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
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
CONFIG="${ROOT}/config/dna_lm_admission.json"
SOURCE="${ROOT}/executions/dna-lm-admission-sources-21064540"
TARGET="${ROOT}/executions/dna-lm-open-safetensors-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite DNA-LM safetensors artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/checkpoints" "${STAGE}/inventory" "${STAGE}/validation"
printf '%s  %s\n' 4e1703b03ff28924904053e842b89079d510a80d2a9ad082de5369adaeea75c8 \
  "${SOURCE}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${SOURCE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

source = Path(sys.argv[1]).resolve(strict=True)
verify_frozen_tree(source)
receipt = json.loads((source / "review" / "source_admission_receipt.json").read_text())
for model_id in ("hyenadna", "caduceus"):
    model = receipt["model_dispositions"][model_id]
    if (
        not model["safe_checkpoint_download_allowed"]
        or model["checkpoint_deserialization_allowed"]
        or model["gpu_forward_allowed"]
    ):
        raise SystemExit(f"safe checkpoint acquisition is not admitted: {model_id}")
print("HyenaDNA/Caduceus source, license, and exposure gates passed")
PY

curl_flags=(--fail --location --retry 8 --retry-all-errors --connect-timeout 30 --continue-at -)
curl "${curl_flags[@]}" \
  --output "${STAGE}/checkpoints/hyenadna.model.safetensors" \
  "https://huggingface.co/LongSafari/hyenadna-large-1m-seqlen-hf/resolve/0a629abf9c7f85b4ec9aa6a1aefa3adcf1907446/model.safetensors"
curl "${curl_flags[@]}" \
  --output "${STAGE}/checkpoints/caduceus.model.safetensors" \
  "https://huggingface.co/kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16/resolve/d89eeb853136ea64da7feb3d0c8e909771b17ae6/model.safetensors"

python "${ROOT}/scripts/safetensors_safe_inventory.py" \
  --input "${STAGE}/checkpoints/hyenadna.model.safetensors" \
  --output "${STAGE}/inventory/hyenadna.json" \
  --expected-sha256 deafb53209bfafb314d14bd81108546a34d473c07ed7ab7a355674376b56ad12 \
  --expected-size 218252304 \
  > "${STAGE}/validation/hyenadna.inventory.stdout.txt" \
  2> "${STAGE}/validation/hyenadna.inventory.stderr.txt"
python "${ROOT}/scripts/safetensors_safe_inventory.py" \
  --input "${STAGE}/checkpoints/caduceus.model.safetensors" \
  --output "${STAGE}/inventory/caduceus.json" \
  --expected-sha256 a3e6976fe90460ff5d90d371b457d0263dc65e156ea5651b4a452e98228daef2 \
  --expected-size 30937760 \
  > "${STAGE}/validation/caduceus.inventory.stdout.txt" \
  2> "${STAGE}/validation/caduceus.inventory.stderr.txt"

python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

root = Path(sys.argv[1]).resolve(strict=True)
models = {}
for model_id in ("hyenadna", "caduceus"):
    inventory = json.loads((root / "inventory" / f"{model_id}.json").read_text())
    if (
        inventory["status"] != "pass"
        or inventory["tensor_data_loaded"]
        or inventory["checkpoint_deserialized"]
        or inventory["tensor_count"] < 1
    ):
        raise SystemExit(f"safe safetensors inventory differs: {model_id}")
    models[model_id] = {
        key: inventory[key]
        for key in (
            "sha256",
            "size_bytes",
            "header_length_bytes",
            "data_length_bytes",
            "tensor_count",
            "parameter_and_buffer_elements",
            "dtype_tensor_counts",
        )
    }
summary = {
    "schema_version": "masld-bench-dna-lm-open-safetensors-admission-v1",
    "status": "pass",
    "models": models,
    "tensor_data_loaded": False,
    "checkpoint_deserialized": False,
    "model_forward_executed": False,
    "observed_outcomes_loaded": False,
    "next_allowed_action": "isolated_no_network_custom_code_runtime_probe_at_common_6000bp_fixture",
}
(root / "inventory" / "admission_summary.json").write_text(
    json.dumps(summary, sort_keys=True, indent=2) + "\n"
)
print(json.dumps(summary, sort_keys=True))
PY

curl --version > "${STAGE}/validation/curl.version.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/safetensors_safe_inventory.py" \
  "${ROOT}/slurm/dna_lm_acquire_open_safetensors.sbatch" \
  "${SOURCE}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads((root / "inventory" / "admission_summary.json").read_text())
freeze_tree(
    root,
    {
        "artifact_class": "dna_lm_open_safetensors",
        "models": ["hyenadna", "caduceus"],
        "checkpoint_sha256": {
            model_id: record["sha256"] for model_id, record in summary["models"].items()
        },
        "tensor_data_loaded": False,
        "checkpoint_deserialized": False,
        "model_forward_executed": False,
        "observed_outcomes_loaded": False,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'DNA-LM open safetensors artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
