#!/usr/bin/env bash
#SBATCH --job-name=model-training
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

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
runtime="${root}/executions/environments/transcriptformer-torch2.5.1-21070738"
runtime_sha256="f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb"
preflight="${root}/executions/transcriptformer-activation-preflight-21070971"
preflight_sha256="2402298e4492a487c063962b4b79d95fe15e871670d972a356a8b24b0a6a588c"
checkpoint="${root}/executions/transcriptformer-materialized-21066178/bundle/tf_sapiens"
python="${runtime}/env/bin/python"
target="${root}/executions/model-training-${SLURM_JOB_ID}"
stage="${target}.staging"
failed="${target}.failed"

preserve_failed_stage() {
  if [[ -d "${stage}" && ! -e "${failed}" ]]; then
    mv "${stage}" "${failed}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${target}" || -e "${stage}" || -e "${failed}" ]]; then
  printf 'Refusing to overwrite model-training artifact.\n' >&2
  exit 1
fi
install -d -m 0750 \
  "${stage}/cache" \
  "${stage}/offline-guard" \
  "${stage}/validation"
printf '%s  %s\n' "${runtime_sha256}" "${runtime}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${preflight_sha256}" "${preflight}/ARTIFACTS.json" \
  | sha256sum --check --strict

cat > "${stage}/offline-guard/sitecustomize.py" <<'PY'
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
export XDG_CACHE_HOME="${stage}/cache"
export MPLCONFIGDIR="${stage}/cache/matplotlib"
export TORCHINDUCTOR_CACHE_DIR="${stage}/cache/torchinductor"
export TRITON_CACHE_DIR="${stage}/cache/triton"
export PYTHONPATH="${stage}/offline-guard:${root}:${root}/src"

"${python}" - <<'PY' > "${stage}/validation/offline-guard.json"
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
  --format=csv,noheader > "${stage}/validation/nvidia-smi.csv"
"${python}" "${root}/scripts/transcriptformer_extract_unlabeled_embeddings.py" \
  --checkpoint-root "${checkpoint}" \
  --fixture-root "${preflight}/fixture" \
  --output "${stage}/embeddings" \
  > "${stage}/validation/embedding.stdout.txt" \
  2> "${stage}/validation/embedding.stderr.txt"
"${python}" -m pip freeze --all > "${stage}/validation/pip-freeze.txt"
sha256sum \
  config/adaptation/transcriptformer.toml \
  scripts/build_transcriptformer_unlabeled_fixture.py \
  scripts/preflight_transcriptformer_cpu.py \
  scripts/transcriptformer_extract_unlabeled_embeddings.py \
  slurm/run_transcriptformer_activation.sbatch \
  "${runtime}/ARTIFACTS.json" \
  "${preflight}/ARTIFACTS.json" \
  > "${stage}/source.sha256"
"${python}" - "${stage}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "embeddings/receipt.json").read_text())
guard = json.loads((root / "validation/offline-guard.json").read_text())
if (
    receipt["status"] != "pass_outcome_blind_embeddings"
    or receipt["rows"] != 1000
    or receipt["embedding_width"] != 2048
    or receipt["evaluation_label_columns_read"]
    or receipt["sealed_outcomes_read"]
    or receipt["downstream_head_fit"]
    or receipt["sealed_champion_eligible"]
    or receipt["exposure_status"] != "unknown"
    or set(receipt["policies"]) != {"native", "common"}
    or not guard["python_INET_socket_blocked"]
):
    raise SystemExit("TranscriptFormer GPU activation firewall differs")
freeze_tree(root, {
    "artifact_class": "transcriptformer_outcome_blind_development_embeddings",
    "dataset_view_id": "resource_atlas_geneformer_smoke_1000_v1",
    "development_rows": 1000,
    "device": "NVIDIA_L40S_sm89",
    "downstream_head_fit": False,
    "evaluation_labels_read": False,
    "exposure_status": "unknown",
    "gene_order_policies": ["native", "common"],
    "sealed_champion_eligible": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${stage}" "${target}"
printf 'model-training artifact: %s\n' "${target}"
sha256sum "${target}/ARTIFACTS.json"
