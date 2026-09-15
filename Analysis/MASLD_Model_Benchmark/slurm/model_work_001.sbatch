#!/usr/bin/env bash
#SBATCH --job-name=model-work-001
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export TOKENIZERS_PARALLELISM=false
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python
CONFIG="${ROOT}/config/gse281364_dna_lm_native_contract.json"
NATIVE="${ROOT}/executions/gse281364-dna-lm-native-contract-21069406"
RUNTIME="${ROOT}/executions/environments/dna-lm-transformers-4.38.1-21064683"
EINOPS="${ROOT}/executions/environments/dnabert2-einops-0.6.1-21064904"
CODE="${ROOT}/executions/dna-lm-custom-code-21064626"
TARGET="${ROOT}/executions/model-work-001-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite model-work-001 artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
declare -a bindings=(
  "92c11e5ae4e71475fff3367d3e5e8ddb801d1d5ea235cfec7a15c1e2251ec623 executions/gse281364-dna-lm-native-contract-21069406"
  "74a3a5c36bd1416daf3fe4b63fe317cfa7e2d008b96e2503d4ca43770b4357de executions/gse281364-dna-lm-common-fixture-21069076"
  "876c072889324661df0500e2f1eb8d09a8f1e3054deddb153d8296b069845a63 executions/environments/dna-lm-transformers-4.38.1-21064683"
  "5b64097862b0db5ad3bdb1d0de1429e932bc06f315e7106209696c702878e3fa executions/dnabert2-runtime-bundle-21064834"
  "25a1d8571cb6e7188f398ac820402ad9f4aea7571a99e5b794cfb92dea57b536 executions/environments/dnabert2-einops-0.6.1-21064904"
  "6c881e7f390d301a818c305bd67e6e305077de87daafd642340814dedbe5f02f executions/dna-lm-weights-only-dnabert2-21064812"
  "4e1703b03ff28924904053e842b89079d510a80d2a9ad082de5369adaeea75c8 executions/dna-lm-admission-sources-21064540"
  "9ec36ffea9bf2c94de1d1d6be6f7f6bc42822a6003f5cf9cd7ee342da2d6f395 executions/dna-lm-weights-only-nucleotide_transformer-21064811"
  "70906faae158909a8972f674f680943d8a51467032ae3266a6a12f2152beb36c executions/dna-lm-custom-code-21064626"
  "6d21df999835624905e08573038dc72e893a5b25a8ca8f4ec1acea8f61d91c77 executions/dna-lm-open-safetensors-21064556"
)
for binding in "${bindings[@]}"; do
  read -r expected relative <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}/ARTIFACTS.json" \
    | sha256sum --check --strict
done
sha256sum --check --strict "${NATIVE}/source.sha256"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${ROOT}" "${bindings[@]}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
for binding in sys.argv[2:]:
    _, relative = binding.split(" ", 1)
    verify_frozen_tree(root / relative)
print("model-work-001 frozen inputs verified")
PY
module purge
install -d -m 0750 \
  "${STAGE}/raw" \
  "${STAGE}/validation" \
  "${STAGE}/hf-home" \
  "${STAGE}/offline-guard"
cat > "${STAGE}/offline-guard/sitecustomize.py" <<'PY'
import socket
_socket = socket.socket
class GuardedSocket(_socket):
    def __new__(cls, family=-1, *args, **kwargs):
        if family in {socket.AF_INET, socket.AF_INET6}:
            raise RuntimeError("network disabled for model-work-001")
        return super().__new__(cls, family, *args, **kwargs)
socket.socket = GuardedSocket
def blocked_connection(*args, **kwargs):
    raise RuntimeError("network disabled for model-work-001")
socket.create_connection = blocked_connection
PY
export HF_HOME="${STAGE}/hf-home"
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export PYTHONPATH="${STAGE}/offline-guard:${ROOT}:${EINOPS}/site-packages:${RUNTIME}/site-packages:${CODE}/code"
"${BASE}" - <<'PY' > "${STAGE}/validation/offline_guard.json"
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
  --format=csv,noheader > "${STAGE}/validation/nvidia-smi.csv"
for model_id in dnabert2 nucleotide_transformer hyenadna; do
  "${BASE}" "${ROOT}/scripts/extract_gse281364_dna_lm_embeddings.py" \
    --root "${ROOT}" \
    --config "${CONFIG}" \
    --model-id "${model_id}" \
    --output "${STAGE}/raw/${model_id}" \
    --seed 20260824 \
    > "${STAGE}/validation/${model_id}.stdout.json" \
    2> "${STAGE}/validation/${model_id}.stderr.txt"
done
"${BASE}" "${ROOT}/scripts/project_gse281364_dna_lm_embeddings.py" \
  --config "${CONFIG}" \
  --raw-root "${STAGE}/raw" \
  --output "${STAGE}/projected" \
  > "${STAGE}/validation/projection.stdout.json" \
  2> "${STAGE}/validation/projection.stderr.txt"
"${BASE}" -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/gse281364_dna_lm_native_contract.py" \
  "${ROOT}/scripts/extract_gse281364_dna_lm_embeddings.py" \
  "${ROOT}/scripts/project_gse281364_dna_lm_embeddings.py" \
  "${ROOT}/slurm/model_work_001.sbatch" \
  "${NATIVE}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
expected = {
    "dnabert2": (768, False, True),
    "nucleotide_transformer": (1280, True, False),
    "hyenadna": (256, False, True),
}
for model, (width, restricted, champion) in expected.items():
    receipt = json.loads((root / f"raw/{model}/receipt.json").read_text())
    if (
        receipt["status"] != "pass_outcome_blind_embedding_extraction"
        or receipt["elements"] != 1033
        or receipt["outer_locus_sequence_groups"] != 1033
        or receipt["outer_folds"] != 5
        or receipt["hidden_width"] != width
        or receipt["restricted_comparator"] != restricted
        or receipt["open_champion_eligible_after_task_gates"] != champion
        or receipt["reporter_counts_read"]
        or receipt["reporter_outcomes_read"]
        or receipt["sealed_labels_read"]
        or receipt["downstream_head_fit"]
        or receipt["runtime_network_allowed"]
    ):
        raise SystemExit(f"{model} embedding receipt differs")
projection = json.loads((root / "projected/receipt.json").read_text())
guard = json.loads((root / "validation/offline_guard.json").read_text())
if (
    projection["status"] != "pass_outcome_blind_fold_safe_projection"
    or projection["elements"] != 1033
    or projection["projection_width"] != 256
    or projection["head_input_width"] != 1024
    or projection["projection_fit_on_held_out_fold"]
    or projection["reporter_counts_read"]
    or projection["reporter_outcomes_read"]
    or projection["sealed_labels_read"]
    or projection["downstream_head_fit"]
    or not guard["python_INET_socket_blocked"]
):
    raise SystemExit("projected head-feature receipt differs")
freeze_tree(root, {
    "artifact_class": "outcome_blind_common_representation_embeddings",
    "dataset_id": "gse281364",
    "models": ["dnabert2", "nucleotide_transformer", "hyenadna"],
    "elements": 1033,
    "outer_locus_sequence_groups": 1033,
    "outer_folds": 5,
    "projection_width": 256,
    "head_input_width": 1024,
    "restricted_comparator": ["nucleotide_transformer"],
    "reporter_outcomes_read": False,
    "downstream_head_fit": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'model-work-001 artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
