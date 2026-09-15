#!/usr/bin/env bash
#SBATCH --job-name=model-training
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
CHECKPOINT=${ROOT}/executions/transcriptformer-materialized-21066178/bundle/tf_sapiens
: "${FIXTURE_ROOT:?FIXTURE_ROOT is required}"
: "${FIXTURE_SHA256:?FIXTURE_SHA256 is required}"
TARGET=${ROOT}/executions/model-training-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite model-training artifact.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/cache" "${STAGE}/offline-guard"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${FIXTURE_SHA256}" "${FIXTURE_ROOT}/ARTIFACTS.json" \
  | sha256sum --check --strict

cat > "${STAGE}/offline-guard/sitecustomize.py" <<'PY'
import socket

_socket = socket.socket
class GuardedSocket(_socket):
    def __new__(cls, family=-1, *args, **kwargs):
        if family in {socket.AF_INET, socket.AF_INET6}:
            raise RuntimeError("network disabled for model-training")
        return super().__new__(cls, family, *args, **kwargs)
socket.socket = GuardedSocket
def blocked_connection(*args, **kwargs):
    raise RuntimeError("network disabled for model-training")
socket.create_connection = blocked_connection
PY
export XDG_CACHE_HOME=${STAGE}/cache
export MPLCONFIGDIR=${STAGE}/cache/matplotlib
export TORCHINDUCTOR_CACHE_DIR=${STAGE}/cache/torchinductor
export TRITON_CACHE_DIR=${STAGE}/cache/triton
export PYTHONPATH=${STAGE}/offline-guard:${ROOT}:${ROOT}/src
cd "${ROOT}"
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader \
  > "${STAGE}/nvidia-smi.csv"
"${RUNTIME}/env/bin/python" -m unittest -v \
  tests.unit.test_transcriptformer_extract_gse256398_shards \
  > "${STAGE}/unittest.txt" 2>&1
"${RUNTIME}/env/bin/python" scripts/transcriptformer_extract_gse256398_shards.py \
  --checkpoint-root "${CHECKPOINT}" \
  --fixture-root "${FIXTURE_ROOT}" \
  --fixture-sha256 "${FIXTURE_SHA256}" \
  --output "${STAGE}/embeddings" \
  > "${STAGE}/embedding.stdout.txt" \
  2> "${STAGE}/embedding.stderr.txt"
"${RUNTIME}/env/bin/python" -m pip freeze --all > "${STAGE}/pip-freeze.txt"
sha256sum \
  scripts/transcriptformer_extract_gse256398_shards.py \
  slurm/run_transcriptformer_gse256398_shards.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  "${FIXTURE_ROOT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
PYTHONPATH=${ROOT}/src "${RUNTIME}/env/bin/python" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "embeddings/embedding_receipt.json").read_text())
if (
    receipt["status"] != "pass_outcome_blind_cross_cohort_embeddings"
    or receipt["rows"] != 172997
    or receipt["donors"] != 26
    or receipt["phenotype_metadata_read"]
    or receipt["barcode_cell_labels_read"]
    or receipt["sealed_outcomes_read"]
    or receipt["supervised_accuracy_claim_allowed"]
    or receipt["sealed_champion_eligible"]
):
    raise SystemExit("GSE256398 TranscriptFormer inference firewall differs")
freeze_tree(root, {
    "artifact_class": "transcriptformer_gse256398_cross_cohort_inference",
    "dataset_id": "gse256398",
    "donors": 26,
    "rows": 172997,
    "device": "NVIDIA_L40S_sm89",
    "phenotype_metadata_read": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
