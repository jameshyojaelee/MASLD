#!/usr/bin/env bash
#SBATCH --job-name=model-check-277
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export MPLBACKEND=Agg
export WANDB_MODE=offline
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome/bin/python
CONFIG="${ROOT}/config/alphagenome_v0_8_sdk_native_task_gate.json"
TARGET="${ROOT}/executions/model-check-277-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
REVISION=71a6beb8c30832f121309a81c2530efa5af7986a
RAW_BASE="https://raw.githubusercontent.com/google-deepmind/alphagenome/${REVISION}"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite AlphaGenome SDK gate: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/sources/github"

curl --fail --silent --show-error --location --proto '=https' --tlsv1.2 \
  https://api.github.com/repos/google-deepmind/alphagenome/releases/latest \
  --output "${STAGE}/sources/github_release_latest.json"
curl --fail --silent --show-error --location --proto '=https' --tlsv1.2 \
  https://api.github.com/repos/google-deepmind/alphagenome/git/ref/tags/v0.8.0 \
  --output "${STAGE}/sources/github_tag_v0.8.0.json"
curl --fail --silent --show-error --location --proto '=https' --tlsv1.2 \
  https://pypi.org/pypi/alphagenome/0.8.0/json \
  --output "${STAGE}/sources/pypi_0.8.0.json"
curl --fail --silent --show-error --location --proto '=https' --tlsv1.2 \
  https://huggingface.co/api/models/google/alphagenome-all-folds \
  --output "${STAGE}/sources/hf_alphagenome_all_folds.json"
curl --fail --silent --show-error --location --proto '=https' --tlsv1.2 \
  "https://huggingface.co/google/alphagenome-all-folds/raw/a8f293a76ee73d5b57f3bf2ae146510589fcf187/README.md" \
  --output "${STAGE}/sources/hf_model_card.md"

while IFS= read -r relative; do
  install -d -m 0750 "${STAGE}/sources/github/$(dirname "${relative}")"
  curl --fail --silent --show-error --location --proto '=https' --tlsv1.2 \
    "${RAW_BASE}/${relative}" \
    --output "${STAGE}/sources/github/${relative}"
done <<'EOF'
README.md
LICENSE
CHANGELOG.md
pyproject.toml
src/alphagenome/models/dna_client.py
src/alphagenome/models/dna_model.py
src/alphagenome/models/dna_output.py
src/alphagenome/models/variant_scorers.py
src/alphagenome/protos/dna_model.proto
src/alphagenome/protos/dna_model_service.proto
EOF

WHEEL_URL=$("${PYTHON}" - "${STAGE}/sources/pypi_0.8.0.json" <<'PY'
import json
from pathlib import Path
import sys

metadata = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
matches = [
    row["url"]
    for row in metadata["urls"]
    if row["filename"] == "alphagenome-0.8.0-py3-none-any.whl"
]
if len(matches) != 1 or not matches[0].startswith("https://files.pythonhosted.org/"):
    raise SystemExit("official wheel URL differs")
print(matches[0])
PY
)
curl --fail --silent --show-error --location --proto '=https' --tlsv1.2 \
  "${WHEEL_URL}" \
  --output "${STAGE}/sources/alphagenome-0.8.0-py3-none-any.whl"

install -d -m 0750 "${STAGE}/validation"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m unittest -v \
  tests.unit.test_audit_alphagenome_v0_8_sdk_native_task_gate \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

"${PYTHON}" "${ROOT}/scripts/audit_alphagenome_v0_8_sdk_native_task_gate.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --source-dir "${STAGE}/sources" \
  --output "${STAGE}/audit" \
  > "${STAGE}/validation/audit.stdout.json" \
  2> "${STAGE}/validation/audit.stderr.txt"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
curl --version > "${STAGE}/validation/curl.version.txt"
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/audit_alphagenome_v0_8_sdk_native_task_gate.py" \
  "${ROOT}/tests/unit/test_audit_alphagenome_v0_8_sdk_native_task_gate.py" \
  "${ROOT}/slurm/audit_alphagenome_v0_8_sdk_native_task_gate_io.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" "${TARGET}" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, publish_directory_noreplace


stage = Path(sys.argv[1]).resolve(strict=True)
target = Path(sys.argv[2])
receipt = json.loads((stage / "audit/receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_sdk_native_interface_admitted_model_execution_blocked"
    or not receipt["sdk_source_admitted"]
    or not receipt["native_interface_admitted"]
    or receipt["offline_runtime_probe"]["sdk_version"] != "0.8.0"
    or receipt["offline_runtime_probe"]["logical_model_version_serialized"] != "ALL_FOLDS"
    or not receipt["offline_runtime_probe"]["grpc_channel_was_fake"]
    or receipt["api_key_inspected"]
    or receipt["api_connection_attempted"]
    or receipt["checkpoint_bytes_downloaded"]
    or receipt["checkpoint_deserialized"]
    or receipt["model_forward_executed"]
    or receipt["observed_outcomes_loaded"]
    or receipt["sealed_assets_loaded"]
    or receipt["metrics_calculated"]
    or receipt["open_champion_eligible"]
    or receipt["conditional_open_model_training_allowed"]
):
    raise SystemExit("AlphaGenome SDK native-task receipt differs")
freeze_tree(
    stage,
    {
        "artifact_class": "restricted_sequence_sdk_native_task_gate",
        "model_id": "alphagenome",
        "sdk_version": "0.8.0",
        "sdk_source_admitted": True,
        "native_interface_admitted": True,
        "api_connection_attempted": False,
        "checkpoint_bytes_downloaded": False,
        "model_forward_executed": False,
        "outcomes_read": False,
        "sealed_assets_read": False,
        "open_champion_eligible": False,
        "status": "pass_sdk_native_interface_admitted_model_execution_blocked",
    },
)
publish_directory_noreplace(stage, target)
PY

trap - EXIT
printf 'AlphaGenome SDK native-task gate: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
