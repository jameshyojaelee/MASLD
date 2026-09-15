#!/usr/bin/env bash
#SBATCH --job-name=model-work-208
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --time=120:00:00
#SBATCH --requeue
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
export PYTHONHASHSEED=20260824
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export OMP_NUM_THREADS=12
export MKL_NUM_THREADS=12

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
CONTRACT=${ROOT}/config/artifacts/models/uce/frozen_screen_50000.json
CONTRACT_SHA256=f771b97b17a8b0be7f8657f1dd3943f0e605844d026ff1101c459edcaa0b4087
: "${UCE_FIXTURE_ROOT:?set UCE_FIXTURE_ROOT to one frozen CPU preflight root}"
: "${UCE_FIXTURE_ARTIFACTS_SHA256:?set exact UCE fixture ARTIFACTS SHA-256}"
: "${MASLD_GPU_DISPATCHER_AUTHORITY:?dispatcher authority is required}"
if [[ "${MASLD_GPU_DISPATCHER_AUTHORITY}" != "sequence_gpu_dispatcher_v2" ]]; then
  printf 'GPU extraction requires the shared dispatcher authority.\n' >&2
  exit 1
fi
MODEL_IDS=${MODEL_IDS:-"uce_4l uce_33l"}
if [[ "${MODEL_IDS}" != "uce_4l uce_33l" ]]; then
  printf 'Frozen UCE bundle must contain 4L then 33L exactly.\n' >&2
  exit 1
fi
UCE_4L_BATCH_SIZE=${UCE_4L_BATCH_SIZE:-32}
UCE_33L_BATCH_SIZE=${UCE_33L_BATCH_SIZE:-8}
ACTIVE_STAGE=
ACTIVE_FAILED=

finish() {
  exit_code=$?
  if [[ "${exit_code}" -ne 0 && -n "${ACTIVE_STAGE}" && -d "${ACTIVE_STAGE}" && ! -e "${ACTIVE_FAILED}" ]]; then
    mv "${ACTIVE_STAGE}" "${ACTIVE_FAILED}"
  fi
  exit "${exit_code}"
}
trap finish EXIT

cd "${ROOT}"
printf '%s  %s\n' \
  "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  "${CONTRACT_SHA256}" "${CONTRACT}" \
  "${UCE_FIXTURE_ARTIFACTS_SHA256}" "${UCE_FIXTURE_ROOT}/ARTIFACTS.json" \
  | sha256sum --check --strict
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${RUNTIME}/env/bin/python" - "${UCE_FIXTURE_ROOT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

for path in map(Path, sys.argv[1:]):
    verify_frozen_tree(path)
PY
nvidia-smi --query-gpu=name,uuid,driver_version,memory.total --format=csv,noheader

declare -A SEEN=()
for MODEL_ID in ${MODEL_IDS}; do
  case "${MODEL_ID}" in
    uce_4l) BATCH_SIZE=${UCE_4L_BATCH_SIZE} ;;
    uce_33l) BATCH_SIZE=${UCE_33L_BATCH_SIZE} ;;
    *) printf 'Unsupported extraction identity: %s\n' "${MODEL_ID}" >&2; exit 1 ;;
  esac
  if [[ -n "${SEEN[${MODEL_ID}]:-}" ]]; then
    printf 'Duplicate extraction identity: %s\n' "${MODEL_ID}" >&2
    exit 1
  fi
  SEEN[${MODEL_ID}]=1
  TARGET=${ROOT}/executions/model-work-208-${MODEL_ID}
  ACTIVE_STAGE=${TARGET}.staging-${SLURM_JOB_ID}
  ACTIVE_FAILED=${TARGET}.failed-${SLURM_JOB_ID}
  if [[ -e "${TARGET}" ]]; then
    "${RUNTIME}/env/bin/python" - "${TARGET}" "${MODEL_ID}" \
      "${UCE_FIXTURE_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import verify_frozen_tree

target = Path(sys.argv[1]).resolve(strict=True)
model_id = sys.argv[2]
fixture_sha256 = sys.argv[3]
manifest = verify_frozen_tree(target)
metadata = manifest.get("metadata", {})
receipt = json.loads((target / "embeddings/receipt.json").read_text(encoding="utf-8"))
if (
    metadata.get("artifact_class") != "uce_frozen_screen_raw_embeddings"
    or metadata.get("model_id") != model_id
    or metadata.get("status") != "passed"
    or metadata.get("downstream_head_fit") is not False
    or receipt.get("model_id") != model_id
    or receipt.get("fixture_artifacts_sha256") != fixture_sha256
    or receipt.get("downstream_head_fit") is not False
    or receipt.get("evaluation_label_columns_read")
    or receipt.get("sealed_outcomes_read") is not False
):
    raise SystemExit("existing UCE model output is not an exact completed chunk")
PY
    printf 'Skipping exact completed immutable model chunk: %s\n' "${MODEL_ID}"
    ACTIVE_STAGE=
    ACTIVE_FAILED=
    continue
  fi
  if [[ -e "${ACTIVE_STAGE}" || -e "${ACTIVE_FAILED}" ]]; then
    printf 'Refusing to overwrite UCE extraction attempt output.\n' >&2
    exit 1
  fi
  "${RUNTIME}/env/bin/python" scripts/extract_uce_frozen_screen_50000.py \
    --project-root "${ROOT}" \
    --contract "${CONTRACT}" \
    --fixture-root "${UCE_FIXTURE_ROOT}" \
    --expected-fixture-artifacts-sha256 "${UCE_FIXTURE_ARTIFACTS_SHA256}" \
    --model-id "${MODEL_ID}" \
    --batch-size "${BATCH_SIZE}" \
    --output "${ACTIVE_STAGE}"
  install -d -m 0750 "${ACTIVE_STAGE}/validation"
  nvidia-smi --query-gpu=name,uuid,driver_version,memory.total \
    --format=csv,noheader > "${ACTIVE_STAGE}/validation/nvidia-smi.csv"
  "${RUNTIME}/env/bin/python" -m pip freeze --all \
    > "${ACTIVE_STAGE}/validation/pip-freeze.txt"
  sha256sum \
    config/artifacts/models/uce/frozen_screen_50000.json \
    scripts/uce_frozen_screen_50000_common.py \
    scripts/extract_uce_frozen_screen_50000.py \
    slurm/run_uce_frozen_screen_50000.sbatch \
    "${RUNTIME}/ARTIFACTS.json" \
    "${UCE_FIXTURE_ROOT}/ARTIFACTS.json" \
    > "${ACTIVE_STAGE}/validation/source.sha256"
  "${RUNTIME}/env/bin/python" - "${ACTIVE_STAGE}" "${MODEL_ID}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
model_id = sys.argv[2]
receipt = json.loads((stage / "embeddings/receipt.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "pass_raw_outcome_blind_extraction"
    or receipt.get("model_id") != model_id
    or receipt.get("rows") != 50_000
    or receipt.get("embedding_width") != 1280
    or "common" not in receipt.get("policies", {})
    or receipt.get("exposure_status") != "target_label_unexposed"
    or receipt.get("evaluation_label_columns_read")
    or receipt.get("sealed_outcomes_read") is not False
    or receipt.get("downstream_head_fit") is not False
):
    raise SystemExit("UCE extraction receipt differs")
freeze_tree(
    stage,
    {
        "artifact_class": "uce_frozen_screen_raw_embeddings",
        "model_id": model_id,
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "development_rows": 50_000,
        "exposure_status": "target_label_unexposed",
        "sealed_champion_eligible": True,
        "evaluation_labels_read": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
        "downstream_head_fit": False,
        "status": "passed",
    },
)
PY
  "${RUNTIME}/env/bin/python" - "${ACTIVE_STAGE}" "${TARGET}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import publish_directory_noreplace

publish_directory_noreplace(Path(sys.argv[1]), Path(sys.argv[2]))
PY
  sha256sum "${TARGET}/ARTIFACTS.json"
  ACTIVE_STAGE=
  ACTIVE_FAILED=
done
trap - EXIT
