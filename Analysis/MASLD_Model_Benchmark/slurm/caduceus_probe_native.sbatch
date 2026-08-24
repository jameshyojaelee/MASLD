#!/usr/bin/env bash
#SBATCH --job-name=masld-caduceus-native
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
RUNTIME="${ROOT}/executions/environments/caduceus-native-torch2.2.0-21065466"
RUNTIME_SHA=3d59792d8a13ca10a2df3d6386bc266eb0e108def6279188793d1cfee4fc1167
CODE="${ROOT}/executions/dna-lm-custom-code-21064626"
CODE_SHA=70906faae158909a8972f674f680943d8a51467032ae3266a6a12f2152beb36c
CHECKPOINTS="${ROOT}/executions/dna-lm-open-safetensors-21064556"
CHECKPOINTS_SHA=6d21df999835624905e08573038dc72e893a5b25a8ca8f4ec1acea8f61d91c77
FIXTURE="${ROOT}/executions/dna-lm-sequence-fixture-21064549"
FIXTURE_SHA=c5018dfd39f4ebde636e76ae4d396bfb3c0b326a8f84894c0b4115f7e8b847b2
TARGET="${ROOT}/executions/caduceus-native-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
RUNTIME_TMP=$(mktemp -d "${TMPDIR:-/tmp}/masld-caduceus-probe-${SLURM_JOB_ID}.XXXXXX")

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed Caduceus native stage: %s\n' "${FAILED}" >&2
  fi
  if [[ -d "${RUNTIME_TMP}" ]]; then
    rm -r "${RUNTIME_TMP}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Caduceus native artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${STAGE}/offline-guard"
for binding in \
  "${RUNTIME_SHA} ${RUNTIME}" \
  "${CODE_SHA} ${CODE}" \
  "${CHECKPOINTS_SHA} ${CHECKPOINTS}" \
  "${FIXTURE_SHA} ${FIXTURE}"; do
  read -r expected directory <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${directory}/ARTIFACTS.json" \
    | sha256sum --check --strict
done

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${RUNTIME}" "${CODE}" "${CHECKPOINTS}" "${FIXTURE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("Caduceus native frozen inputs verified")
PY
module purge

python3 - "${RUNTIME}/validation/dependency_metadata.json" \
  "${STAGE}/validation/dependency_terms_review.json" <<'PY'
import json
from pathlib import Path
import sys

records = json.loads(Path(sys.argv[1]).read_text())
expected = {
    "mamba-ssm": "1.2.0.post1",
    "causal-conv1d": "1.2.0.post2",
    "torch": "2.2.0",
    "transformers": "4.38.1",
    "safetensors": "0.4.2",
}
observed = {record["name"]: record for record in records}
if any(
    name not in observed
    or observed[name]["version"] != version
    or not any(
        marker in " ".join(
            [str(observed[name]["license_metadata"])]
            + observed[name]["license_classifiers"]
        ).lower()
        for marker in ("apache", "bsd", "mit", "mpl")
    )
    for name, version in expected.items()
):
    raise SystemExit("Caduceus runtime dependency terms metadata differs")
receipt = {
    "status": "pass",
    "review_scope": "runtime_dependency_metadata_before_checkpoint_execution",
    "required_dependencies": expected,
    "model_code_license": "Apache-2.0",
    "checkpoint_card_license": "apache-2.0",
    "derivative_checkpoint_redistribution": "not_exercised_by_probe",
}
Path(sys.argv[2]).write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
print("Caduceus runtime dependency terms metadata review passed")
PY

cat > "${STAGE}/offline-guard/sitecustomize.py" <<'PY'
import socket

_socket = socket.socket

class GuardedSocket(_socket):
    def __new__(cls, family=-1, *args, **kwargs):
        if family in {socket.AF_INET, socket.AF_INET6}:
            raise RuntimeError("network disabled for Caduceus checkpoint execution")
        return super().__new__(cls, family, *args, **kwargs)

socket.socket = GuardedSocket

def blocked_connection(*args, **kwargs):
    raise RuntimeError("network disabled for Caduceus checkpoint execution")

socket.create_connection = blocked_connection
PY

tar -xzf "${RUNTIME}/payload/caduceus-native-env.tar.gz" -C "${RUNTIME_TMP}"
export HF_HOME="${STAGE}/hf-home"
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export PYTHONPATH="${STAGE}/offline-guard"
install -d -m 0750 "${HF_HOME}"
"${RUNTIME_TMP}/bin/python" "${ROOT}/scripts/dna_lm_probe_caduceus.py" \
  --code-root "${CODE}/code/caduceus" \
  --checkpoint "${CHECKPOINTS}/checkpoints/caduceus.model.safetensors" \
  --fixture "${FIXTURE}" \
  --output "${STAGE}/predictions" \
  > "${STAGE}/validation/probe.stdout.txt" \
  2> "${STAGE}/validation/probe.stderr.txt"
"${RUNTIME_TMP}/bin/python" -m pip freeze --all \
  > "${STAGE}/validation/effective-pip-freeze.txt"
nvidia-smi -q > "${STAGE}/validation/nvidia-smi-q.txt"
sha256sum \
  "${ROOT}/scripts/dna_lm_probe_caduceus.py" \
  "${ROOT}/slurm/caduceus_probe_native.sbatch" \
  "${RUNTIME}/ARTIFACTS.json" \
  "${CODE}/ARTIFACTS.json" \
  "${CHECKPOINTS}/ARTIFACTS.json" \
  "${FIXTURE}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
rm -r "${RUNTIME_TMP}"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "predictions/runtime_probe_receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["learned_parameter_count"] != 7_725_312
    or receipt["checkpoint_parameter_and_buffer_elements"] != 7_725_344
    or receipt["checkpoint_tensor_count"] != 276
    or receipt["fixture_count"] != 3
    or receipt["feature_width"] != 256
    or not receipt["tied_alias_parameter_identity_validated"]
    or receipt["deterministic_masked_repeat_abs_difference"] > 1e-5
    or receipt["deterministic_feature_repeat_max_abs_difference"] > 1e-5
    or not receipt["input_gradient_backward_finite"]
    or receipt["input_gradient_max_abs"] <= 0.0
    or not receipt["compiled_mamba_forward_executed"]
    or receipt["head_fit"]
    or receipt["observed_outcomes_loaded"]
    or receipt["sealed_outcomes_loaded"]
    or receipt["runtime_network_allowed"]
):
    raise SystemExit("Caduceus native probe receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "caduceus_native_runtime_probe",
        "runtime_artifacts_sha256": "3d59792d8a13ca10a2df3d6386bc266eb0e108def6279188793d1cfee4fc1167",
        "code_artifacts_sha256": "70906faae158909a8972f674f680943d8a51467032ae3266a6a12f2152beb36c",
        "checkpoint_artifacts_sha256": "6d21df999835624905e08573038dc72e893a5b25a8ca8f4ec1acea8f61d91c77",
        "fixture_artifacts_sha256": "c5018dfd39f4ebde636e76ae4d396bfb3c0b326a8f84894c0b4115f7e8b847b2",
        "checkpoint_sha256": "a3e6976fe90460ff5d90d371b457d0263dc65e156ea5651b4a452e98228daef2",
        "terms": "Apache-2.0_open_comparator_and_task_candidate",
        "outcome_firewall": "no_observed_or_sealed_outcomes",
        "head_fit": False,
        "champion_eligible_only_after_task_and_external_evaluation_gates": True,
    },
)
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Caduceus native runtime probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
