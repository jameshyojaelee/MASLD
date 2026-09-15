#!/usr/bin/env bash
#SBATCH --job-name=model-probe-043
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
export TOKENIZERS_PARALLELISM=false
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/geneformer-v2-316m-compat-runtime-21077474
RUNTIME_ARTIFACTS_SHA256=c034ee71fb2c67f04f11e91b41ea8c9bf6a7b9df3a4bb72378e5b2696058c307
RUNTIME_ARCHIVE_SHA256=a4fdc1f1d0b67c4753354c5a04ccc05c84a3bf7442d8655bc8ed6493a2487775
BUNDLE=${ROOT}/executions/geneformer-acquisition-21014827/bundles/geneformer_v2_316m--11f311a5078fc4c8
INVENTORY=${ROOT}/executions/geneformer_v2_316m-safe-inventory-21066211/checkpoint_inventory.json
TOKENS=${ROOT}/executions/geneformer-sweep/21015001/geneformer_v2_316m__two_layer_mlp/prepare/rank_value_tokens.jsonl
TOKEN_SUMMARY=${ROOT}/executions/geneformer-sweep/21015001/geneformer_v2_316m__two_layer_mlp/prepare/tokenization_summary.json
PRIOR_EMBEDDINGS=${ROOT}/executions/embeddings/geneformer_v2_316m/embeddings.npy
PRIOR_ROW_ORDER=${ROOT}/executions/embeddings/geneformer_v2_316m/embedding_row_order.txt
ADAPTER=${ROOT}/src/masld_bench/adapters/geneformer_common_lane.py
EXTRACTOR=${ROOT}/executions/fiveseed/extract_geneformer_emb.py
TARGET=${ROOT}/executions/geneformer-v2-316m-one-batch-l40s-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
WORK=$(mktemp -d "${TMPDIR:-/tmp}/geneformer-l40s.${SLURM_JOB_ID}.XXXXXX")

finish() {
  exit_code=$?
  if [[ "${exit_code}" -ne 0 && -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed one-batch preflight: %s\n' "${FAILED}" >&2
  fi
  if [[ -d "${WORK}" ]]; then
    rm -rf -- "${WORK}"
  fi
  exit "${exit_code}"
}
trap finish EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite one-batch preflight: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/offline-guard" "${STAGE}/validation" "${WORK}/runtime"
cd "${ROOT}"

printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${RUNTIME_ARCHIVE_SHA256}" \
  "${RUNTIME}/payload/geneformer-decima-compat-env.tar.gz" \
  | sha256sum --check --strict

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${RUNTIME}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

manifest = verify_frozen_tree(Path(sys.argv[1]))
metadata = manifest.get("metadata", {})
if (
    metadata.get("status") != "passed_cpu_gate"
    or metadata.get("strict_checkpoint_restore") is not True
    or metadata.get("raw_unrefined_embeddings_only") is not True
    or metadata.get("prior_smoke_numeric_reproduction_completed") is not False
):
    raise SystemExit("Geneformer compatibility runtime gate differs")
PY
module purge

tar -xzf "${RUNTIME}/payload/geneformer-decima-compat-env.tar.gz" -C "${WORK}/runtime"
cat > "${STAGE}/offline-guard/sitecustomize.py" <<'PY'
import socket

_socket = socket.socket


class GuardedSocket(_socket):
    def __init__(self, family=-1, *args, **kwargs):
        if family in (socket.AF_INET, socket.AF_INET6):
            raise RuntimeError("network disabled by Geneformer preflight")
        super().__init__(family, *args, **kwargs)


socket.socket = GuardedSocket
PY
export PYTHONPATH="${STAGE}/offline-guard"
"${WORK}/runtime/bin/python" - <<'PY'
import socket

try:
    socket.socket(socket.AF_INET, socket.SOCK_STREAM)
except RuntimeError as error:
    if str(error) != "network disabled by Geneformer preflight":
        raise
else:
    raise SystemExit("offline network guard did not activate")
PY

nvidia-smi --query-gpu=name,uuid,driver_version,memory.total \
  --format=csv,noheader > "${STAGE}/validation/nvidia-smi.csv"
"${WORK}/runtime/bin/python" \
  "${ROOT}/scripts/preflight_geneformer_v2_316m_common.py" \
  --mode gpu \
  --bundle "${BUNDLE}" \
  --inventory "${INVENTORY}" \
  --tokens "${TOKENS}" \
  --token-summary "${TOKEN_SUMMARY}" \
  --prior-embeddings "${PRIOR_EMBEDDINGS}" \
  --prior-row-order "${PRIOR_ROW_ORDER}" \
  --adapter "${ADAPTER}" \
  --extraction-script "${EXTRACTOR}" \
  --output "${STAGE}/gpu_gate_receipt.json" \
  --output-embeddings "${STAGE}/raw_embeddings.npy"

sha256sum \
  "${RUNTIME}/ARTIFACTS.json" \
  "${RUNTIME}/payload/geneformer-decima-compat-env.tar.gz" \
  "${ROOT}/scripts/preflight_geneformer_v2_316m_common.py" \
  "${ROOT}/slurm/preflight_geneformer_v2_316m_one_batch_l40s.sbatch" \
  > "${STAGE}/validation/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${RUNTIME_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
runtime_artifacts_sha256 = sys.argv[2]
receipt = json.loads((stage / "gpu_gate_receipt.json").read_text(encoding="utf-8"))
result = receipt.get("result", {})
if (
    receipt.get("status") != "pass_l40s_prior_batch_reproduction"
    or receipt.get("raw_unrefined_embeddings_only") is not True
    or receipt.get("labels_loaded") is not False
    or receipt.get("sealed_outcomes_loaded") is not False
    or result.get("prior_smoke_raw_embeddings_reproduced") is not True
):
    raise SystemExit("Geneformer one-batch reproduction receipt differs")
freeze_tree(
    stage,
    {
        "artifact_class": "geneformer_v2_316m_l40s_one_batch_preflight",
        "model_id": "geneformer_v2_316m",
        "runtime_artifacts_sha256": runtime_artifacts_sha256,
        "raw_unrefined_embeddings_only": True,
        "batch_rows": result["batch_rows"],
        "prior_max_abs": result["prior_max_abs"],
        "prior_mean_abs": result["prior_mean_abs"],
        "repeat_max_abs": result["repeat_max_abs"],
        "prior_smoke_raw_embeddings_reproduced": True,
        "project_labels_read": False,
        "sealed_outcomes_read": False,
        "status": "passed_l40s_preflight",
    },
)
PY

trap - EXIT
rm -rf -- "${WORK}"
mv "${STAGE}" "${TARGET}"
printf 'Geneformer V2-316M L40S preflight: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
