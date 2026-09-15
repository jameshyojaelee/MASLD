#!/usr/bin/env bash
#SBATCH --job-name=model-training-500
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8
: "${MASLD_GPU_DISPATCHER_AUTHORITY:?dispatcher authority is required}"
if [[ "${MASLD_GPU_DISPATCHER_AUTHORITY}" != sequence_gpu_dispatcher_v2 ]]; then
  printf 'Invalid GPU dispatcher authority.\n' >&2
  exit 2
fi

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
CHECKPOINT=${ROOT}/executions/transcriptformer-materialized-21066178/bundle/tf_sapiens
FIXTURE=${ROOT}/executions/model-data-051-21077382
FIXTURE_ARTIFACTS_SHA256=ba40ce9555bb6d7bc278be2a80015b993f9f3650fa06f59505ad279c88cc117a
PREFLIGHT=${ROOT}/executions/model-check-032-21077481
PREFLIGHT_ARTIFACTS_SHA256=de8dd839b1bbb7d5110370a56b8484c47320c08360cfd24ea3a1f8c3c5f68308
PYTHON=${RUNTIME}/env/bin/python
TARGET=${ROOT}/executions/model-training-500-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite model-training artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/cache" \
  "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${FIXTURE_ARTIFACTS_SHA256}" "${FIXTURE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${PREFLIGHT_ARTIFACTS_SHA256}" "${PREFLIGHT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' \
  cdcc828a24c7176e08b002a59bb80fde1fe89e29d7ba6fadc00cc58a1e3b1a48 \
  "${ROOT}/scripts/transcriptformer_extract_unlabeled_embeddings.py" \
  | sha256sum --check --strict
printf '%s  %s\n' \
  a726ae976c1f1c23a64a4da4e3084a3192fac3522d023ac6f920ab81dff00c9a \
  "${ROOT}/tests/unit/test_transcriptformer_extract_unlabeled_embeddings.py" \
  | sha256sum --check --strict
printf '%s  %s\n' \
  ed4bd726db8210d11c034b806c54fa541afce62664ee0846912022036e414635 \
  "${ROOT}/scripts/offline_guard/sitecustomize.py" \
  | sha256sum --check --strict

export XDG_CACHE_HOME="${STAGE}/cache"
export MPLCONFIGDIR="${STAGE}/cache/matplotlib"
export TORCHINDUCTOR_CACHE_DIR="${STAGE}/cache/torchinductor"
export TRITON_CACHE_DIR="${STAGE}/cache/triton"
export PYTHONPATH="${ROOT}/scripts/offline_guard:${ROOT}:${ROOT}/src"

cd "${ROOT}"
"${PYTHON}" -m unittest \
  tests.unit.test_transcriptformer_extract_unlabeled_embeddings \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${PYTHON}" - "${FIXTURE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

value = verify_frozen_tree(Path(sys.argv[1]))
if (
    value["metadata"].get("artifact_class")
    != "transcriptformer_frozen_screen_input"
    or value["metadata"].get("development_rows") != 50_000
    or value["metadata"].get("evaluation_labels_read")
    or value["metadata"].get("sealed_outcomes_read")
):
    raise SystemExit("frozen-screen input firewall differs")
PY
"${PYTHON}" - <<'PY' > "${STAGE}/validation/offline_guard.json"
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
nvidia-smi --query-gpu=name,driver_version,memory.total \
  --format=csv,noheader > "${STAGE}/validation/nvidia_smi.csv"
"${PYTHON}" scripts/transcriptformer_extract_unlabeled_embeddings.py \
  --checkpoint-root "${CHECKPOINT}" \
  --fixture-root "${FIXTURE}/fixture" \
  --output "${STAGE}/embeddings" \
  --expected-rows 50000 \
  > "${STAGE}/validation/embedding_stdout.txt" \
  2> "${STAGE}/validation/embedding_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/transcriptformer_extract_unlabeled_embeddings.py \
  scripts/offline_guard/sitecustomize.py \
  tests/unit/test_transcriptformer_extract_unlabeled_embeddings.py \
  slurm/run_model_training_500.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  "${FIXTURE}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "embeddings/receipt.json").read_text(encoding="utf-8"))
guard = json.loads((root / "validation/offline_guard.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "pass_outcome_blind_embeddings"
    or receipt.get("rows") != 50_000
    or receipt.get("embedding_width") != 2_048
    or set(receipt.get("policies", {})) != {"native", "common"}
    or receipt.get("evaluation_label_columns_read")
    or receipt.get("sealed_outcomes_read")
    or receipt.get("downstream_head_fit")
    or receipt.get("exposure_status") != "unknown"
    or not guard.get("python_INET_socket_blocked")
):
    raise SystemExit("50k model embedding firewall differs")
freeze_tree(root, {
    "artifact_class": "cell_foundation_frozen_screen_embeddings",
    "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
    "development_rows": 50_000,
    "downstream_head_fit": False,
    "evaluation_labels_read": False,
    "exposure_status": "unknown",
    "gene_order_policies": ["native", "common"],
    "model_id": "transcriptformer_tf_sapiens",
    "sealed_champion_eligible": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
