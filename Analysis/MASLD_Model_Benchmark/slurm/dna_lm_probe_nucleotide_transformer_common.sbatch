#!/usr/bin/env bash
#SBATCH --job-name=masld-nt-common
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
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
SOURCE="${ROOT}/executions/dna-lm-admission-sources-21064540"
SOURCE_SHA=4e1703b03ff28924904053e842b89079d510a80d2a9ad082de5369adaeea75c8
WEIGHTS="${ROOT}/executions/dna-lm-weights-only-nucleotide_transformer-21064811"
WEIGHTS_SHA=9ec36ffea9bf2c94de1d1d6be6f7f6bc42822a6003f5cf9cd7ee342da2d6f395
FIXTURE="${ROOT}/executions/dna-lm-sequence-fixture-21064549"
FIXTURE_SHA=c5018dfd39f4ebde636e76ae4d396bfb3c0b326a8f84894c0b4115f7e8b847b2
TARGET="${ROOT}/executions/nucleotide-transformer-common-runtime-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed NT stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite NT probe artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/validation" \
  "${STAGE}/hf-home" \
  "${STAGE}/offline-guard"
for binding in \
  "${RUNTIME_SHA} ${RUNTIME}" \
  "${SOURCE_SHA} ${SOURCE}" \
  "${WEIGHTS_SHA} ${WEIGHTS}" \
  "${FIXTURE_SHA} ${FIXTURE}"; do
  read -r expected directory <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${directory}/ARTIFACTS.json" \
    | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${RUNTIME}" "${SOURCE}" "${WEIGHTS}" "${FIXTURE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("Nucleotide Transformer frozen inputs verified")
PY
module purge

cat > "${STAGE}/offline-guard/sitecustomize.py" <<'PY'
import socket

_socket = socket.socket

class GuardedSocket(_socket):
    def __new__(cls, family=-1, *args, **kwargs):
        if family in {socket.AF_INET, socket.AF_INET6}:
            raise RuntimeError("network disabled for checkpoint execution")
        return super().__new__(cls, family, *args, **kwargs)

socket.socket = GuardedSocket

def blocked_connection(*args, **kwargs):
    raise RuntimeError("network disabled for checkpoint execution")

socket.create_connection = blocked_connection
PY
export HF_HOME="${STAGE}/hf-home"
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export TOKENIZERS_PARALLELISM=false
export PYTHONPATH="${STAGE}/offline-guard:${RUNTIME}/site-packages"
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
"${BASE}" "${ROOT}/scripts/dna_lm_probe_nucleotide_transformer.py" \
  --config "${SOURCE}/sources/nt_hf_config" \
  --vocab "${SOURCE}/sources/nt_hf_vocab" \
  --tokenizer-config "${SOURCE}/sources/nt_hf_tokenizer_config" \
  --checkpoint "${WEIGHTS}/converted/nucleotide_transformer.model.safetensors" \
  --fixture "${FIXTURE}" \
  --output "${STAGE}/predictions" \
  --seed 20260824 \
  > "${STAGE}/validation/probe.stdout.txt" \
  2> "${STAGE}/validation/probe.stderr.txt"
"${BASE}" -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/dna_lm_probe_nucleotide_transformer.py" \
  "${ROOT}/slurm/dna_lm_probe_nucleotide_transformer_common.sbatch" \
  "${RUNTIME}/ARTIFACTS.json" \
  "${SOURCE}/ARTIFACTS.json" \
  "${WEIGHTS}/ARTIFACTS.json" \
  "${FIXTURE}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "predictions" / "runtime_probe_receipt.json").read_text())
guard = json.loads((root / "validation" / "offline_guard.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["head_fit"]
    or receipt["observed_outcomes_loaded"]
    or receipt["sealed_outcomes_loaded"]
    or receipt["open_champion_eligible"]
    or not guard["python_INET_socket_blocked"]
):
    raise SystemExit("NT runtime receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "nucleotide_transformer_common_runtime_probe",
        "runtime_artifacts_sha256": "876c072889324661df0500e2f1eb8d09a8f1e3054deddb153d8296b069845a63",
        "source_artifacts_sha256": "4e1703b03ff28924904053e842b89079d510a80d2a9ad082de5369adaeea75c8",
        "weights_artifacts_sha256": "9ec36ffea9bf2c94de1d1d6be6f7f6bc42822a6003f5cf9cd7ee342da2d6f395",
        "fixture_artifacts_sha256": "c5018dfd39f4ebde636e76ae4d396bfb3c0b326a8f84894c0b4115f7e8b847b2",
        "terms": "cc-by-nc-sa-4.0_restricted_comparator",
        "outcome_firewall": "development_sequence_only_no_observed_or_sealed_outcomes",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'NT common runtime probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
