#!/usr/bin/env bash
#SBATCH --job-name=model-work-302
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

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
runtime_root="${root}/executions/environments/scprint2-native-torch2.8.0-v2-21073597"
runtime="${runtime_root}/source/.venv"
fixture="${root}/executions/scprint2-no-neighbor-fixture-21075645"
preflight="${root}/executions/scprint2-native-runtime-preflight-v2-21074027"
checkpoint="${root}/executions/scprint2-acquisition-21066044/18hebyht-final-small.ckpt"
target="${root}/executions/model-work-302-${SLURM_JOB_ID}"
stage="${target}.staging"
failed="${target}.failed"

preserve_failed_stage() {
  if [[ -d "${stage}" && ! -e "${failed}" ]]; then
    mv "${stage}" "${failed}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${target}" || -e "${stage}" || -e "${failed}" ]]; then
  printf 'Refusing to overwrite model-work-302 artifact.\n' >&2
  exit 1
fi
install -d -m 0750 "${stage}/cache" "${stage}/offline-guard" "${stage}/validation"

printf '%s  %s\n' \
  "4474a0d069d09a9a805e7d6e33936af95eeab9fb72d6e3c131e7e7733350842d" \
  "${runtime_root}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' \
  "edf71c212340ebfdcc266ef7c0ecb73449193b36a17147e53d3801804082c9f7" \
  "${fixture}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' \
  "9ed6d6db96bca280b2b2a99f8c3256be329978ad145e234aaac4815ef33295fa" \
  "${preflight}/ARTIFACTS.json" | sha256sum --check --strict

cat > "${stage}/offline-guard/sitecustomize.py" <<'PY'
import socket

_socket = socket.socket
class GuardedSocket(_socket):
    def __new__(cls, family=-1, *args, **kwargs):
        if family in {socket.AF_INET, socket.AF_INET6}:
            raise RuntimeError("network disabled for model-work-302")
        return super().__new__(cls, family, *args, **kwargs)
socket.socket = GuardedSocket
def blocked_connection(*args, **kwargs):
    raise RuntimeError("network disabled for model-work-302")
socket.create_connection = blocked_connection
PY

export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export XDG_CACHE_HOME="${stage}/cache"
export MPLCONFIGDIR="${stage}/cache/matplotlib"
export TORCHINDUCTOR_CACHE_DIR="${stage}/cache/torchinductor"
export TRITON_CACHE_DIR="${stage}/cache/triton"
export LD_LIBRARY_PATH="/gpfs/commons/home/jameslee/micromamba/envs/scprint2x/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export PYTHONPATH="${stage}/offline-guard:${root}:${root}/src"

"${runtime}/bin/python" - <<'PY' > "${stage}/validation/offline-guard.json"
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
"${runtime}/bin/python" "${root}/scripts/scprint2_extract_no_neighbor_embeddings.py" \
  --checkpoint "${checkpoint}" \
  --fixture-root "${fixture}" \
  --runtime-root "${runtime_root}" \
  --preflight-root "${preflight}" \
  --output "${stage}/embeddings" \
  > "${stage}/validation/embedding.stdout.txt" \
  2> "${stage}/validation/embedding.stderr.txt"
"${runtime}/bin/python" - <<'PY' > "${stage}/validation/pip-freeze.txt"
from importlib.metadata import distributions

records = {
    f"{dist.metadata.get('Name', 'UNKNOWN')}=={dist.version}"
    for dist in distributions()
}
print("\n".join(sorted(records, key=str.casefold)))
PY
sha256sum \
  scripts/scprint2_extract_no_neighbor_embeddings.py \
  slurm/model_work_302.sbatch \
  "${runtime_root}/ARTIFACTS.json" \
  "${fixture}/ARTIFACTS.json" \
  "${preflight}/ARTIFACTS.json" \
  > "${stage}/source.sha256"

"${runtime}/bin/python" - "${stage}" <<'PY'
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
    or receipt["donors"] != 102
    or receipt["neighbor_tensors_used"]
    or receipt["evaluation_label_columns_read"]
    or receipt["sealed_outcomes_read"]
    or receipt["sealed_champion_eligible"]
    or not receipt["source_prefilter_depth_used"]
    or receipt["exposure_status"] != "unknown"
    or set(receipt["policies"]) != {
        "common", "random_checkpoint_max", "released_selection_default_no_neighbor"
    }
    or not guard["python_INET_socket_blocked"]
):
    raise SystemExit("scPRINT-2 outcome-blind embedding firewall differs")
freeze_tree(root, {
    "artifact_class": "scprint2_outcome_blind_development_embeddings",
    "dataset_view_id": "resource_atlas_geneformer_smoke_1000_v1",
    "development_rows": 1000,
    "device": "NVIDIA_L40S_sm89",
    "evaluation_labels_read": False,
    "exposure_status": "unknown",
    "model_id": "scprint2_small_v2",
    "neighbor_tensors_used": False,
    "sealed_champion_eligible": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY

trap - EXIT
mv "${stage}" "${target}"
sha256sum "${target}/ARTIFACTS.json"
