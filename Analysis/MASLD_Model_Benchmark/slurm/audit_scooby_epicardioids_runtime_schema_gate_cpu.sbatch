#!/usr/bin/env bash
#SBATCH --job-name=model-check-279
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
export OMP_NUM_THREADS=2
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/borzoi/bin/python
CONFIG="${ROOT}/config/scooby_epicardioids_runtime_schema_gate.json"
TARGET="${ROOT}/executions/model-check-279-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
REVISION=dbfe6cd150d018d46a6d24ffa0406d7f2e89b471
RAW_BASE="https://raw.githubusercontent.com/gagneurlab/scooby/${REVISION}"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Scooby schema gate: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/sources/scooby/modeling" "${STAGE}/validation"

curl --fail --silent --show-error --location --proto '=https' --tlsv1.2 \
  "${RAW_BASE}/LICENSE" --output "${STAGE}/sources/LICENSE"
curl --fail --silent --show-error --location --proto '=https' --tlsv1.2 \
  "${RAW_BASE}/setup.py" --output "${STAGE}/sources/setup.py"
curl --fail --silent --show-error --location --proto '=https' --tlsv1.2 \
  "${RAW_BASE}/scooby/modeling/scooby.py" --output "${STAGE}/sources/scooby/modeling/scooby.py"
curl --fail --silent --show-error --location --proto '=https' --tlsv1.2 \
  https://huggingface.co/api/models/lauradmartens/epicardioids-scooby \
  --output "${STAGE}/sources/hf_api.json"

export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m unittest -v \
  tests.unit.test_audit_scooby_epicardioids_runtime_schema_gate \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${PYTHON}" "${ROOT}/scripts/audit_scooby_epicardioids_runtime_schema_gate.py" \
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
  "${ROOT}/scripts/audit_scooby_epicardioids_runtime_schema_gate.py" \
  "${ROOT}/tests/unit/test_audit_scooby_epicardioids_runtime_schema_gate.py" \
  "${ROOT}/slurm/audit_scooby_epicardioids_runtime_schema_gate_cpu.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" "${TARGET}" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, publish_directory_noreplace


stage = Path(sys.argv[1]).resolve(strict=True)
target = Path(sys.argv[2])
receipt = json.loads((stage / "audit/receipt.json").read_text(encoding="utf-8"))
probe = receipt["runtime_probe"]
if (
    receipt["status"] != "pass_checkpoint_download_and_strict_restore_probe_allowed"
    or probe["state_dict_tensors"] != 225
    or probe["state_dict_elements_including_buffers"] != 193241454
    or probe["schema_missing"]
    or probe["schema_unexpected"]
    or probe["schema_mismatched"]
    or not probe["deterministic_repeat_bit_identical"]
    or probe["checkpoint_weights_loaded"]
    or probe["full_sequence_forward_executed"]
    or receipt["checkpoint_payload_downloaded"]
    or receipt["checkpoint_deserialized"]
    or receipt["checkpoint_forward_executed"]
    or receipt["project_data_read"]
    or receipt["sealed_data_read"]
    or receipt["rna_conditioned_atac_eligible"]
    or receipt["open_champion_eligible"]
):
    raise SystemExit("Scooby runtime schema receipt differs")
freeze_tree(
    stage,
    {
        "artifact_class": "restricted_observed_atac_model_runtime_schema_gate",
        "model_id": "scooby_epicardioids",
        "checkpoint_payload_downloaded": False,
        "checkpoint_forward_executed": False,
        "project_data_read": False,
        "sealed_data_read": False,
        "rna_conditioned_atac_eligible": False,
        "open_champion_eligible": False,
        "status": "pass_checkpoint_download_and_strict_restore_probe_allowed",
    },
)
publish_directory_noreplace(stage, target)
PY

trap - EXIT
printf 'Scooby runtime schema gate: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
