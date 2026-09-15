#!/usr/bin/env bash
#SBATCH --job-name=model-training-503
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
: "${MASLD_GPU_DISPATCHER_AUTHORITY:?dispatcher authority is required}"
: "${ONE_BATCH_ROOT:?passed one-batch artifact root is required}"
: "${ONE_BATCH_ARTIFACTS_SHA256:?passed one-batch artifact hash is required}"
: "${SCREEN_PREFLIGHT_ROOT:?50k CPU preflight artifact root is required}"
: "${SCREEN_PREFLIGHT_ARTIFACTS_SHA256:?50k CPU preflight artifact hash is required}"
if [[ "${MASLD_GPU_DISPATCHER_AUTHORITY}" != sequence_gpu_dispatcher_v2 ]]; then
  printf 'Invalid GPU dispatcher authority.\n' >&2
  exit 1
fi
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
export TOKENIZERS_PARALLELISM=false
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export OMP_NUM_THREADS=6
export MKL_NUM_THREADS=6

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/geneformer-v2-316m-compat-runtime-21077474
RUNTIME_SHA256=c034ee71fb2c67f04f11e91b41ea8c9bf6a7b9df3a4bb72378e5b2696058c307
RUNTIME_ARCHIVE_SHA256=a4fdc1f1d0b67c4753354c5a04ccc05c84a3bf7442d8655bc8ed6493a2487775
FIXTURE=${ROOT}/executions/geneformer-v2-316m-frozen-screen-fixture-21078339
FIXTURE_SHA256=bb94045a46c7231d50a17ffd921a146ad6858b39bcfcda516da61230fcfb9bf6
BUNDLE=${ROOT}/executions/geneformer-acquisition-21014827/bundles/geneformer_v2_316m--11f311a5078fc4c8
BUNDLE_SHA256=b269cef53fa92cc7e7ab9c21a2b85de69904cdfaaf5881e9676c57dc2ae02ae0
TARGET=${ROOT}/executions/geneformer-v2-316m-frozen-screen-extraction-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
WORK=$(mktemp -d "${TMPDIR:-/tmp}/geneformer-50k-extract.${SLURM_JOB_ID}.XXXXXX")

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
trap finish EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Geneformer 50k extraction.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/offline-guard" "${STAGE}/validation" "${WORK}/runtime"
cd "${ROOT}"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${ROOT}" "${RUNTIME}" "${RUNTIME_SHA256}" \
  "${FIXTURE}" "${FIXTURE_SHA256}" \
  "${BUNDLE}" "${BUNDLE_SHA256}" \
  "${ONE_BATCH_ROOT}" "${ONE_BATCH_ARTIFACTS_SHA256}" \
  "${SCREEN_PREFLIGHT_ROOT}" "${SCREEN_PREFLIGHT_ARTIFACTS_SHA256}" <<'PY'
from hashlib import sha256
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
pairs = [(Path(path), digest) for path, digest in zip(sys.argv[2::2], sys.argv[3::2], strict=True)]
manifests = []
for path, expected in pairs:
    resolved = path.resolve(strict=True)
    try:
        resolved.relative_to(root)
    except ValueError as error:
        raise SystemExit("required artifact escapes benchmark root") from error
    observed = sha256((resolved / "ARTIFACTS.json").read_bytes()).hexdigest()
    if observed != expected:
        raise SystemExit("required artifact hash differs")
    manifests.append(verify_frozen_tree(resolved))
one_batch = manifests[3].get("metadata", {})
screen = manifests[4].get("metadata", {})
if (
    one_batch.get("status") != "passed_l40s_preflight"
    or one_batch.get("prior_smoke_raw_embeddings_reproduced") is not True
    or one_batch.get("project_labels_read") is not False
    or screen.get("status") != "passed"
    or screen.get("full_extraction_submitted") is not False
    or screen.get("evaluation_labels_read") is not False
):
    raise SystemExit("GPU activation gates differ")
PY
module purge
printf '%s  %s\n' "${RUNTIME_ARCHIVE_SHA256}" \
  "${RUNTIME}/payload/geneformer-decima-compat-env.tar.gz" \
  | sha256sum --check --strict
tar -xzf "${RUNTIME}/payload/geneformer-decima-compat-env.tar.gz" -C "${WORK}/runtime"

cat > "${STAGE}/offline-guard/sitecustomize.py" <<'PY'
import socket

_socket = socket.socket


class GuardedSocket(_socket):
    def __init__(self, family=-1, *args, **kwargs):
        if family in (socket.AF_INET, socket.AF_INET6):
            raise RuntimeError("network disabled by Geneformer extraction")
        super().__init__(family, *args, **kwargs)


socket.socket = GuardedSocket
PY
export PYTHONPATH="${STAGE}/offline-guard:${ROOT}"
nvidia-smi --query-gpu=name,uuid,driver_version,memory.total \
  --format=csv,noheader > "${STAGE}/validation/nvidia-smi-before.csv"
"${WORK}/runtime/bin/python" \
  scripts/extract_geneformer_v2_316m_frozen_screen_50000.py \
  --fixture-root "${FIXTURE}/fixture" \
  --model-directory "${BUNDLE}/upstream/Geneformer-V2-316M" \
  --output "${STAGE}/extraction" \
  --expected-rows 50000 \
  --batch-size 32 \
  > "${STAGE}/validation/extraction.stdout.txt" \
  2> "${STAGE}/validation/extraction.stderr.txt"
nvidia-smi --query-gpu=name,uuid,driver_version,memory.total \
  --format=csv,noheader > "${STAGE}/validation/nvidia-smi-after.csv"
sha256sum \
  scripts/preflight_geneformer_v2_316m_common.py \
  scripts/extract_geneformer_v2_316m_frozen_screen_50000.py \
  slurm/run_geneformer_v2_316m_frozen_screen_50000.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  "${FIXTURE}/ARTIFACTS.json" \
  "${BUNDLE}/ARTIFACTS.json" \
  "${ONE_BATCH_ROOT}/ARTIFACTS.json" \
  "${SCREEN_PREFLIGHT_ROOT}/ARTIFACTS.json" \
  > "${STAGE}/validation/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" \
  "${ONE_BATCH_ARTIFACTS_SHA256}" "${SCREEN_PREFLIGHT_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((stage / "extraction/extraction_receipt.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "pass_raw_outcome_blind_extraction"
    or receipt.get("rows") != 50_000
    or receipt.get("embedding_shape") != [50_000, 1152]
    or receipt.get("raw_unrefined_embeddings_only") is not True
    or receipt.get("downstream_head_fit") is not False
    or receipt.get("evaluation_labels_read") is not False
    or receipt.get("histology_read") is not False
    or receipt.get("sealed_outcomes_read") is not False
):
    raise SystemExit("Geneformer 50k extraction receipt differs")
freeze_tree(
    stage,
    {
        "artifact_class": "geneformer_v2_316m_frozen_screen_raw_embeddings",
        "model_id": "geneformer_v2_316m",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "rows": 50_000,
        "one_batch_artifacts_sha256": sys.argv[2],
        "screen_preflight_artifacts_sha256": sys.argv[3],
        "raw_unrefined_embeddings_only": True,
        "downstream_head_fit": False,
        "evaluation_labels_read": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
PY

trap - EXIT
rm -rf -- "${WORK}"
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
