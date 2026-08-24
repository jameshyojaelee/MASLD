#!/usr/bin/env bash
#SBATCH --job-name=masld-dna-lm-pickle
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=48:00:00
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
TARGET="${ROOT}/executions/dna-lm-safe-pickle-checkpoints-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite DNA-LM pickle-checkpoint artifact: %s\n' "${TARGET}" >&2
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
for model_id in ("dnabert2", "nucleotide_transformer", "evo"):
    model = receipt["model_dispositions"][model_id]
    if (
        not model["safe_checkpoint_download_allowed"]
        or model["checkpoint_deserialization_allowed"]
        or model["gpu_forward_allowed"]
    ):
        raise SystemExit(f"safe checkpoint acquisition is not admitted: {model_id}")
print("DNABERT-2/NT/Evo source, license, and exposure gates passed")
PY

curl_flags=(--fail --location --retry 8 --retry-all-errors --connect-timeout 30 --continue-at -)
curl "${curl_flags[@]}" \
  --output "${STAGE}/checkpoints/dnabert2.pytorch_model.bin" \
  "https://huggingface.co/zhihan1996/DNABERT-2-117M/resolve/7bce263b15377fc15361f52cfab88f8b586abda0/pytorch_model.bin"
curl "${curl_flags[@]}" \
  --output "${STAGE}/checkpoints/nucleotide_transformer.pytorch_model.bin" \
  "https://huggingface.co/InstaDeepAI/nucleotide-transformer-500m-human-ref/resolve/f87b5d7233295242e79c951873d290f4cf992045/pytorch_model.bin"
curl "${curl_flags[@]}" \
  --output "${STAGE}/checkpoints/evo.evo2_7b.pt" \
  "https://huggingface.co/arcinstitute/evo2_7b/resolve/bda0089f92582d5baabf0f22d9fc85f3588f6b58/evo2_7b.pt"

python "${ROOT}/scripts/torch_checkpoint_safe_inventory.py" \
  --input "${STAGE}/checkpoints/dnabert2.pytorch_model.bin" \
  --output "${STAGE}/inventory/dnabert2.json" \
  --expected-sha256 7ff39ec77a484dd01070a41bfd6e95cdd7247bec80fe357ab43a4be33687aeba \
  --expected-size 468354983 \
  > "${STAGE}/validation/dnabert2.inventory.stdout.txt" \
  2> "${STAGE}/validation/dnabert2.inventory.stderr.txt"
python "${ROOT}/scripts/torch_checkpoint_safe_inventory.py" \
  --input "${STAGE}/checkpoints/nucleotide_transformer.pytorch_model.bin" \
  --output "${STAGE}/inventory/nucleotide_transformer.json" \
  --expected-sha256 971cb721bd8cc8134d3665b5d94b195fffb18a1480ad5925d6dcc9290bf1c17f \
  --expected-size 1942934337 \
  > "${STAGE}/validation/nucleotide_transformer.inventory.stdout.txt" \
  2> "${STAGE}/validation/nucleotide_transformer.inventory.stderr.txt"
python "${ROOT}/scripts/torch_checkpoint_safe_inventory.py" \
  --input "${STAGE}/checkpoints/evo.evo2_7b.pt" \
  --output "${STAGE}/inventory/evo.json" \
  --expected-sha256 c66645929dc1b9c631f5be656da8726f38946315dc9167000a615dd626fcecf4 \
  --expected-size 13766621200 \
  --maximum-uncompressed-bytes 25769803776 \
  > "${STAGE}/validation/evo.inventory.stdout.txt" \
  2> "${STAGE}/validation/evo.inventory.stderr.txt"

python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

root = Path(sys.argv[1]).resolve(strict=True)
models = {}
for model_id in ("dnabert2", "nucleotide_transformer", "evo"):
    inventory = json.loads((root / "inventory" / f"{model_id}.json").read_text())
    if (
        inventory["status"] != "pass"
        or inventory["archive_extracted"]
        or inventory["pickle_imported"]
        or inventory["torch_imported"]
        or inventory["checkpoint_deserialized"]
        or inventory["tensor_storage_member_count"] < 1
    ):
        raise SystemExit(f"safe torch ZIP inventory differs: {model_id}")
    models[model_id] = {
        key: inventory[key]
        for key in (
            "sha256",
            "size_bytes",
            "container_type",
            "member_count",
            "regular_file_count",
            "uncompressed_regular_bytes",
            "pickle_metadata_members",
            "tensor_storage_member_count",
        )
    }
summary = {
    "schema_version": "masld-bench-dna-lm-pickle-checkpoint-admission-v1",
    "status": "pass",
    "models": models,
    "archive_extracted": False,
    "pickle_imported": False,
    "torch_imported": False,
    "checkpoint_deserialized": False,
    "model_forward_executed": False,
    "observed_outcomes_loaded": False,
    "next_allowed_action": "isolated_no_network_weights_only_state_dict_conversion_with_key_shape_dtype_freeze",
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
  "${ROOT}/scripts/torch_checkpoint_safe_inventory.py" \
  "${ROOT}/slurm/dna_lm_acquire_pickle_checkpoints.sbatch" \
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
        "artifact_class": "dna_lm_safe_pickle_checkpoints",
        "models": ["dnabert2", "nucleotide_transformer", "evo"],
        "checkpoint_sha256": {
            model_id: record["sha256"] for model_id, record in summary["models"].items()
        },
        "archive_extracted": False,
        "pickle_imported": False,
        "torch_imported": False,
        "checkpoint_deserialized": False,
        "model_forward_executed": False,
        "observed_outcomes_loaded": False,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'DNA-LM safe pickle-checkpoint artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
