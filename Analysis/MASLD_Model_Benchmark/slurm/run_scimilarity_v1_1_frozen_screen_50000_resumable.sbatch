#!/usr/bin/env bash
#SBATCH --job-name=model-training-504
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --requeue
#SBATCH --signal=B:USR1@120
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
: "${MASLD_GPU_DISPATCHER_AUTHORITY:?dispatcher authority is required}"
if [[ "${MASLD_GPU_DISPATCHER_AUTHORITY}" != "sequence_gpu_dispatcher_v2" ]]; then
  printf 'Invalid GPU dispatcher authority.\n' >&2
  exit 1
fi
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/scimilarity-v1-1-compat-runtime-21079070
RUNTIME_SHA256=8f66e2e478742665072a23ed4b115346830fcee516094c59f20164c94b3a5a4a
RUNTIME_ARCHIVE_SHA256=e2d0a97c8631d3cdacb4234fa7bc167ddd8fc391c3dd2a2b7dd23ebaebe88548
FIXTURE=${ROOT}/executions/scimilarity-v1-1-frozen-screen-fixture-21079292
FIXTURE_SHA256=b6d50557107b5e6a878480741ec6a4c781663218876943ae205578b1a4ecd1e5
PREFLIGHT=${ROOT}/executions/scimilarity-v1-1-frozen-screen-preflight-21079643
PREFLIGHT_SHA256=8dcf5c8d91f59ebfa5cbb9d55be7d3299c0eeaae6522b9be7c3a893306eeb6e9
QUEUE_VALIDATION=${ROOT}/executions/model-training-504-queue-validation-v1
MODEL=${ROOT}/executions/scimilarity-model-20260823/model_v1.1
STATE=${ROOT}/executions/resumable-model-training-504/scimilarity-v1-1-frozen-screen-v1
TARGET=${ROOT}/executions/scimilarity-v1-1-frozen-screen-extraction-v1
RESTART_COUNT=${SLURM_RESTART_COUNT:-0}
STAGE=${TARGET}.staging-${SLURM_JOB_ID}-${RESTART_COUNT}
FAILED=${TARGET}.failed-${SLURM_JOB_ID}-${RESTART_COUNT}
PAUSED=${TARGET}.paused-${SLURM_JOB_ID}-${RESTART_COUNT}
WORK=$(mktemp -d "${TMPDIR:-/tmp}/model-training-504.${SLURM_JOB_ID}.${RESTART_COUNT}.XXXXXX")
CHILD_PID=

finish() {
  exit_code=$?
  if [[ "${exit_code}" -ne 0 && -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
  if [[ -d "${WORK}" ]]; then
    rm -rf -- "${WORK}"
  fi
  exit "${exit_code}"
}
forward_stop() {
  if [[ -n "${CHILD_PID}" ]] && kill -0 "${CHILD_PID}" 2>/dev/null; then
    kill -USR1 "${CHILD_PID}"
  fi
}
trap finish EXIT
trap forward_stop USR1

cd "${ROOT}"
module purge
module load python/3.11.5-GCCcore-13.2.0
CONTROL_PYTHON=$(command -v python)

if [[ -d "${TARGET}" ]]; then
  PYTHONPATH="${ROOT}/src" "${CONTROL_PYTHON}" - "${TARGET}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

target = Path(sys.argv[1]).resolve(strict=True)
manifest = verify_frozen_tree(target)
receipt = json.loads((target / "extraction/receipt.json").read_text(encoding="utf-8"))
if (
    manifest.get("metadata", {}).get("status") != "passed"
    or receipt.get("status") != "pass_raw_outcome_blind_resumable_embeddings"
    or receipt.get("rows") != 50_000
    or receipt.get("downstream_head_fit") is not False
    or receipt.get("metrics_run") is not False
    or receipt.get("sealed_outcomes_read") is not False
):
    raise SystemExit("completed SCimilarity extraction differs")
PY
  trap - EXIT
  rm -rf -- "${WORK}"
  printf 'Completed frozen extraction already exists: %s\n' "${TARGET}"
  exit 0
fi
if [[ -e "${STAGE}" || -e "${FAILED}" || -e "${PAUSED}" ]]; then
  printf 'SCimilarity extraction attempt path already exists.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/offline-guard" "${STAGE}/validation" "${WORK}/runtime"
install -d -m 0750 "${STATE}"

PYTHONPATH="${ROOT}/src" "${CONTROL_PYTHON}" - \
  "${RUNTIME}" "${RUNTIME_SHA256}" \
  "${FIXTURE}" "${FIXTURE_SHA256}" \
  "${PREFLIGHT}" "${PREFLIGHT_SHA256}" \
  "${QUEUE_VALIDATION}" <<'PY'
from hashlib import sha256
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

values = sys.argv[1:]
for path_value, expected in zip(values[:6:2], values[1:6:2], strict=True):
    path = Path(path_value).resolve(strict=True)
    observed = sha256((path / "ARTIFACTS.json").read_bytes()).hexdigest()
    if observed != expected:
        raise SystemExit("required SCimilarity artifact hash differs")
    verify_frozen_tree(path)
validation = verify_frozen_tree(Path(values[6]).resolve(strict=True))
if validation.get("metadata", {}).get("bundle_id") != "model-training-504":
    raise SystemExit("SCimilarity queue validation differs")
PY
printf '%s  %s\n' "${RUNTIME_ARCHIVE_SHA256}" \
  "${RUNTIME}/payload/scimilarity-v1-1-compat-env.tar.gz" \
  | sha256sum --check --strict
tar -xzf "${RUNTIME}/payload/scimilarity-v1-1-compat-env.tar.gz" -C "${WORK}/runtime"

cat > "${STAGE}/offline-guard/sitecustomize.py" <<'PY'
import socket

_socket = socket.socket


class GuardedSocket(_socket):
    def __init__(self, family=-1, *args, **kwargs):
        if family in (socket.AF_INET, socket.AF_INET6):
            raise RuntimeError("network disabled by frozen embedding extraction")
        super().__init__(family, *args, **kwargs)


socket.socket = GuardedSocket
PY
export PYTHONPATH="${STAGE}/offline-guard:${ROOT}:${ROOT}/src"
nvidia-smi --query-gpu=name,uuid,driver_version,memory.total \
  --format=csv,noheader > "${STAGE}/validation/nvidia-smi-before.csv"

set +e
"${WORK}/runtime/bin/python" \
  scripts/extract_scimilarity_v1_1_frozen_screen_50000_resumable.py \
  --fixture-root "${FIXTURE}/fixture" \
  --checkpoint "${MODEL}/encoder.ckpt" \
  --state-root "${STATE}" \
  --output "${STAGE}/extraction" \
  --expected-rows 50000 \
  --expected-dimension 28231 \
  --batch-size 512 \
  --chunk-size 4096 \
  --device cuda \
  > "${STAGE}/validation/extraction.stdout.txt" \
  2> "${STAGE}/validation/extraction.stderr.txt" &
CHILD_PID=$!
wait_status=0
while kill -0 "${CHILD_PID}" 2>/dev/null; do
  wait "${CHILD_PID}" || wait_status=$?
done
if wait "${CHILD_PID}" 2>/dev/null; then
  wait_status=0
elif [[ "${wait_status}" -eq 0 ]]; then
  wait_status=$?
fi
CHILD_PID=
set -e

if [[ "${wait_status}" -eq 75 ]]; then
  mv "${STAGE}" "${PAUSED}"
  trap - EXIT
  rm -rf -- "${WORK}"
  scontrol requeue "${SLURM_JOB_ID}"
  exit 0
fi
if [[ "${wait_status}" -ne 0 ]]; then
  exit "${wait_status}"
fi

nvidia-smi --query-gpu=name,uuid,driver_version,memory.total \
  --format=csv,noheader > "${STAGE}/validation/nvidia-smi-after.csv"
"${WORK}/runtime/bin/python" -m pip freeze --all \
  > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  scripts/extract_scimilarity_v1_1_frozen_screen_50000.py \
  scripts/extract_scimilarity_v1_1_frozen_screen_50000_resumable.py \
  slurm/run_scimilarity_v1_1_frozen_screen_50000_resumable.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  "${FIXTURE}/ARTIFACTS.json" \
  "${PREFLIGHT}/ARTIFACTS.json" \
  "${QUEUE_VALIDATION}/ARTIFACTS.json" \
  > "${STAGE}/validation/source.sha256"

PYTHONPATH="${ROOT}/src" "${CONTROL_PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((stage / "extraction/receipt.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "pass_raw_outcome_blind_resumable_embeddings"
    or receipt.get("rows") != 50_000
    or receipt.get("embedding_shape") != [50_000, 128]
    or receipt.get("aggregate_development_exposure") != "encoder_seen"
    or receipt.get("study_gene_observability_shape") != [7, 28_231]
    or receipt.get("raw_unrefined_embeddings_only") is not True
    or receipt.get("released_reference_index_used") is not False
    or receipt.get("downstream_head_fit") is not False
    or receipt.get("common_head_fit") is not False
    or receipt.get("metrics_run") is not False
    or receipt.get("evaluation_label_columns_read")
    or receipt.get("histology_read") is not False
    or receipt.get("sealed_outcomes_read") is not False
    or receipt.get("sealed_champion_eligible") is not False
):
    raise SystemExit("SCimilarity full extraction receipt differs")
freeze_tree(
    stage,
    {
        "artifact_class": "scimilarity_v1_1_frozen_screen_raw_embeddings",
        "model_id": "scimilarity_v1_1",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "development_rows": 50_000,
        "aggregate_development_exposure": "encoder_seen",
        "resumable_chunk_state_used": True,
        "raw_unrefined_embeddings_only": True,
        "released_reference_index_used": False,
        "downstream_head_fit": False,
        "common_head_fit": False,
        "metrics_run": False,
        "evaluation_labels_read": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
        "sealed_champion_eligible": False,
        "status": "passed",
    },
)
PY
PYTHONPATH="${ROOT}/src" "${CONTROL_PYTHON}" - "${STAGE}" "${TARGET}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import publish_directory_noreplace

publish_directory_noreplace(Path(sys.argv[1]), Path(sys.argv[2]))
PY

trap - EXIT
rm -rf -- "${WORK}"
sha256sum "${TARGET}/ARTIFACTS.json"
