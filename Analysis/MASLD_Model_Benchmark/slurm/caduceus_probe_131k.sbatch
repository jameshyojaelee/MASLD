#!/usr/bin/env bash
#SBATCH --job-name=masld-caduceus-131k
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/caduceus-native-torch2.2.0-21065466"
RUNTIME_SHA=3d59792d8a13ca10a2df3d6386bc266eb0e108def6279188793d1cfee4fc1167
CODE="${ROOT}/executions/dna-lm-custom-code-21064626"
CODE_SHA=70906faae158909a8972f674f680943d8a51467032ae3266a6a12f2152beb36c
CHECKPOINTS="${ROOT}/executions/dna-lm-open-safetensors-21064556"
CHECKPOINTS_SHA=6d21df999835624905e08573038dc72e893a5b25a8ca8f4ec1acea8f61d91c77
FIXTURE="${ROOT}/executions/dna-lm-sequence-fixture-21064549"
FIXTURE_SHA=c5018dfd39f4ebde636e76ae4d396bfb3c0b326a8f84894c0b4115f7e8b847b2
TARGET="${ROOT}/executions/caduceus-131k-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
RUNTIME_TMP=$(mktemp -d "${TMPDIR:-/tmp}/masld-caduceus-131k-${SLURM_JOB_ID}.XXXXXX")

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed Caduceus 131-kb stage: %s\n' "${FAILED}" >&2
  fi
  if [[ -d "${RUNTIME_TMP}" ]]; then rm -r "${RUNTIME_TMP}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Caduceus 131-kb artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${STAGE}/offline-guard"
for binding in \
  "${RUNTIME_SHA} ${RUNTIME}" "${CODE_SHA} ${CODE}" \
  "${CHECKPOINTS_SHA} ${CHECKPOINTS}" "${FIXTURE_SHA} ${FIXTURE}"; do
  read -r expected directory <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${directory}/ARTIFACTS.json" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${RUNTIME}" "${CODE}" "${CHECKPOINTS}" "${FIXTURE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]: verify_frozen_tree(Path(value))
print("Caduceus 131-kb frozen inputs verified")
PY
module purge

cat > "${STAGE}/offline-guard/sitecustomize.py" <<'PY'
import socket
_socket = socket.socket
class GuardedSocket(_socket):
    def __new__(cls, family=-1, *args, **kwargs):
        if family in {socket.AF_INET, socket.AF_INET6}:
            raise RuntimeError("network disabled for Caduceus 131-kb execution")
        return super().__new__(cls, family, *args, **kwargs)
socket.socket = GuardedSocket
socket.create_connection = lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("network disabled"))
PY
tar -xzf "${RUNTIME}/payload/caduceus-native-env.tar.gz" -C "${RUNTIME_TMP}"
export HF_HOME="${STAGE}/hf-home"
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export PYTHONPATH="${STAGE}/offline-guard"
install -d -m 0750 "${HF_HOME}"
"${RUNTIME_TMP}/bin/python" "${ROOT}/scripts/dna_lm_probe_caduceus_131k.py" \
  --code-root "${CODE}/code/caduceus" \
  --checkpoint "${CHECKPOINTS}/checkpoints/caduceus.model.safetensors" \
  --fixture "${FIXTURE}" \
  --manifest "${FIXTURE}/fixture/sequence_manifest.tsv" \
  --fasta "${FIXTURE}/fixture/caduceus_native_131072.alleles.fa.gz" \
  --output "${STAGE}/predictions" \
  > "${STAGE}/validation/probe.stdout.txt" \
  2> "${STAGE}/validation/probe.stderr.txt"
"${RUNTIME_TMP}/bin/python" -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
nvidia-smi -q > "${STAGE}/validation/nvidia-smi-q.txt"
sha256sum \
  "${ROOT}/scripts/dna_lm_probe_caduceus_131k.py" \
  "${ROOT}/slurm/caduceus_probe_131k.sbatch" \
  "${RUNTIME}/ARTIFACTS.json" "${CODE}/ARTIFACTS.json" \
  "${CHECKPOINTS}/ARTIFACTS.json" "${FIXTURE}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
rm -r "${RUNTIME_TMP}"

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
    or receipt["input_length_bp"] != 131_072
    or receipt["profile_count"] != 4
    or receipt["feature_width"] != 256
    or receipt["numeric_repeat_abs_difference"] > 1e-5
    or receipt["head_fit"]
    or receipt["observed_outcomes_loaded"]
    or receipt["sealed_outcomes_loaded"]
    or receipt["runtime_network_allowed"]
): raise SystemExit("Caduceus 131-kb receipt differs")
freeze_tree(root, {
    "artifact_class": "caduceus_native_131k_runtime_probe",
    "runtime_artifacts_sha256": "3d59792d8a13ca10a2df3d6386bc266eb0e108def6279188793d1cfee4fc1167",
    "checkpoint_sha256": "a3e6976fe90460ff5d90d371b457d0263dc65e156ea5651b4a452e98228daef2",
    "fixture_artifacts_sha256": "c5018dfd39f4ebde636e76ae4d396bfb3c0b326a8f84894c0b4115f7e8b847b2",
    "outcome_firewall": "no_observed_or_sealed_outcomes",
    "head_fit": False,
})
PY
module purge
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Caduceus 131-kb native probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
