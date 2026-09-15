#!/usr/bin/env bash
#SBATCH --job-name=model-data-281
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/borzoi/bin/python
CONFIG="${ROOT}/config/scooby_epicardioids_checkpoint_acquisition.json"
FAILED_INPUT="${ROOT}/executions/model-data-280-21097262.failed"
TARGET="${ROOT}/executions/model-data-281-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Scooby checkpoint acquisition v2: %s\n' "${TARGET}" >&2
  exit 1
fi
if [[ ! -f "${FAILED_INPUT}/checkpoint/model.safetensors" ]]; then
  printf 'Failed-attempt payload is absent: %s\n' "${FAILED_INPUT}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/checkpoint" "${STAGE}/validation"

if [[ ! -f "${FAILED_INPUT}/COMPLETE" ]]; then
  export PYTHONPATH="${ROOT}:${ROOT}/src"
  "${PYTHON}" - "${FAILED_INPUT}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree


failed = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    failed,
    {
        "artifact_class": "failed_restricted_model_checkpoint_acquisition",
        "source_job_id": "21097262",
        "failure": "local_bounded_reader_overread_by_one_byte",
        "checkpoint_payload_sha256": "f945e1937a46ddd5cb31249ef71e721e7a7b1381aa207fb2de449f23ac75a9e8",
        "checkpoint_payload_size_bytes": 772988628,
        "checkpoint_deserialized": False,
        "model_forward_executed": False,
        "project_data_read": False,
        "sealed_data_read": False,
        "status": "failed_preserved",
    },
)
PY
fi
sha256sum "${FAILED_INPUT}/ARTIFACTS.json" > "${STAGE}/validation/failed_attempt_artifacts.sha256"
printf '%s  %s\n' \
  "f945e1937a46ddd5cb31249ef71e721e7a7b1381aa207fb2de449f23ac75a9e8" \
  "${FAILED_INPUT}/checkpoint/model.safetensors" | sha256sum --check --strict
cp --reflink=auto "${FAILED_INPUT}/checkpoint/model.safetensors" "${STAGE}/checkpoint/model.safetensors"

export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m unittest -v \
  tests.unit.test_acquire_scooby_epicardioids_checkpoint \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${PYTHON}" "${ROOT}/scripts/acquire_scooby_epicardioids_checkpoint.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --checkpoint "${STAGE}/checkpoint/model.safetensors" \
  --output "${STAGE}/audit" \
  > "${STAGE}/validation/audit.stdout.json" \
  2> "${STAGE}/validation/audit.stderr.txt"

"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
cp "${FAILED_INPUT}/ARTIFACTS.json" "${STAGE}/validation/failed_attempt_ARTIFACTS.json"
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/acquire_scooby_epicardioids_checkpoint.py" \
  "${ROOT}/tests/unit/test_acquire_scooby_epicardioids_checkpoint.py" \
  "${ROOT}/slurm/acquire_scooby_epicardioids_checkpoint_io_v2.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" "${TARGET}" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, publish_directory_noreplace


stage = Path(sys.argv[1]).resolve(strict=True)
target = Path(sys.argv[2])
receipt = json.loads((stage / "audit/receipt.json").read_text(encoding="utf-8"))
checkpoint = receipt["checkpoint"]
if (
    receipt["status"] != "pass_exact_checkpoint_acquired_safe_restore_allowed"
    or checkpoint["sha256"] != "f945e1937a46ddd5cb31249ef71e721e7a7b1381aa207fb2de449f23ac75a9e8"
    or checkpoint["size_bytes"] != 772988628
    or checkpoint["state_dict_tensors"] != 225
    or checkpoint["state_dict_elements_including_buffers"] != 193241454
    or not checkpoint["full_payload_hash_verified"]
    or not receipt["checkpoint_payload_downloaded"]
    or receipt["checkpoint_deserialized"]
    or receipt["model_forward_executed"]
    or receipt["full_sequence_forward_executed"]
    or receipt["project_data_read"]
    or receipt["sealed_data_read"]
    or receipt["task_scope"]["rna_conditioned_atac_eligible"]
    or receipt["license"]["open_champion_eligible"]
):
    raise SystemExit("Scooby checkpoint acquisition v2 receipt differs")
freeze_tree(
    stage,
    {
        "artifact_class": "restricted_observed_atac_model_exact_checkpoint",
        "model_id": "scooby_epicardioids",
        "source_failed_job_id": "21097262",
        "checkpoint_sha256": checkpoint["sha256"],
        "checkpoint_payload_downloaded": True,
        "checkpoint_deserialized": False,
        "model_forward_executed": False,
        "project_data_read": False,
        "sealed_data_read": False,
        "rna_conditioned_atac_eligible": False,
        "open_champion_eligible": False,
        "status": "pass_exact_checkpoint_acquired_safe_restore_allowed",
    },
)
publish_directory_noreplace(stage, target)
PY

trap - EXIT
printf 'Scooby exact checkpoint acquisition v2: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${FAILED_INPUT}/ARTIFACTS.json"
