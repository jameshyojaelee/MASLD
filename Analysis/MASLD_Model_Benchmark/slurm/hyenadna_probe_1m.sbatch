#!/usr/bin/env bash
#SBATCH --job-name=masld-hyenadna-1m
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python
RUNTIME="${ROOT}/executions/environments/dna-lm-transformers-4.38.1-21064683"
RUNTIME_SHA=876c072889324661df0500e2f1eb8d09a8f1e3054deddb153d8296b069845a63
CODE="${ROOT}/executions/dna-lm-custom-code-21064626"
CODE_SHA=70906faae158909a8972f674f680943d8a51467032ae3266a6a12f2152beb36c
WEIGHTS="${ROOT}/executions/dna-lm-open-safetensors-21064556"
WEIGHTS_SHA=6d21df999835624905e08573038dc72e893a5b25a8ca8f4ec1acea8f61d91c77
FIXTURE="${ROOT}/executions/dna-lm-sequence-fixture-21064549"
FIXTURE_SHA=c5018dfd39f4ebde636e76ae4d396bfb3c0b326a8f84894c0b4115f7e8b847b2
TARGET="${ROOT}/executions/hyenadna-1m-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed HyenaDNA 1-Mb stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite HyenaDNA 1-Mb artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${STAGE}/hf-home" "${STAGE}/offline-guard"
for binding in \
  "${RUNTIME_SHA} ${RUNTIME}" "${CODE_SHA} ${CODE}" \
  "${WEIGHTS_SHA} ${WEIGHTS}" "${FIXTURE_SHA} ${FIXTURE}"; do
  read -r expected directory <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${directory}/ARTIFACTS.json" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${RUNTIME}" "${CODE}" "${WEIGHTS}" "${FIXTURE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]: verify_frozen_tree(Path(value))
print("HyenaDNA 1-Mb frozen inputs verified")
PY
module purge
cat > "${STAGE}/offline-guard/sitecustomize.py" <<'PY'
import socket
_socket = socket.socket
class GuardedSocket(_socket):
    def __new__(cls, family=-1, *args, **kwargs):
        if family in {socket.AF_INET, socket.AF_INET6}:
            raise RuntimeError("network disabled for HyenaDNA 1-Mb execution")
        return super().__new__(cls, family, *args, **kwargs)
socket.socket = GuardedSocket
def blocked_connection(*args, **kwargs): raise RuntimeError("network disabled")
socket.create_connection = blocked_connection
PY
export HF_HOME="${STAGE}/hf-home"
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export TOKENIZERS_PARALLELISM=false
export PYTHONPATH="${STAGE}/offline-guard:${RUNTIME}/site-packages:${CODE}/code"
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader \
  > "${STAGE}/validation/nvidia-smi.csv"
"${BASE}" "${ROOT}/scripts/dna_lm_probe_hyenadna_1m.py" \
  --code-root "${CODE}/code/hyenadna" \
  --checkpoint "${WEIGHTS}/checkpoints/hyenadna.model.safetensors" \
  --fixture "${FIXTURE}" \
  --manifest "${FIXTURE}/fixture/sequence_manifest.tsv" \
  --fasta "${FIXTURE}/fixture/hyenadna_native_1000000.alleles.fa.gz" \
  --output "${STAGE}/predictions" \
  > "${STAGE}/validation/probe.stdout.txt" \
  2> "${STAGE}/validation/probe.stderr.txt"
"${BASE}" -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/dna_lm_probe_hyenadna_1m.py" \
  "${ROOT}/slurm/hyenadna_probe_1m.sbatch" \
  "${RUNTIME}/ARTIFACTS.json" "${CODE}/ARTIFACTS.json" \
  "${WEIGHTS}/ARTIFACTS.json" "${FIXTURE}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "predictions/probe_receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["input_length_bp"] != 1_000_000
    or receipt["forward_count"] != 5
    or receipt["numeric_score_repeat_abs_difference"] > 1e-6
    or receipt["numeric_feature_repeat_max_abs_difference"] > 1e-5
    or receipt["head_fit"]
    or receipt["observed_outcomes_loaded"]
    or receipt["sealed_outcomes_loaded"]
    or receipt["runtime_network_allowed"]
): raise SystemExit("HyenaDNA 1-Mb receipt differs")
freeze_tree(root, {
    "artifact_class": "hyenadna_native_1m_runtime_probe",
    "runtime_artifacts_sha256": "876c072889324661df0500e2f1eb8d09a8f1e3054deddb153d8296b069845a63",
    "weights_sha256": "deafb53209bfafb314d14bd81108546a34d473c07ed7ab7a355674376b56ad12",
    "fixture_artifacts_sha256": "c5018dfd39f4ebde636e76ae4d396bfb3c0b326a8f84894c0b4115f7e8b847b2",
    "outcome_firewall": "no_observed_or_sealed_outcomes",
    "head_fit": False,
})
PY
module purge
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'HyenaDNA 1-Mb native probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
