#!/usr/bin/env bash
#SBATCH --job-name=model-check-282
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/borzoi/bin/python
CONFIG="${ROOT}/config/scooby_epicardioids_checkpoint_decoder_probe.json"
TARGET="${ROOT}/executions/model-check-282-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Scooby checkpoint decoder probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m unittest -v \
  tests.unit.test_probe_scooby_epicardioids_checkpoint_decoder \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${PYTHON}" "${ROOT}/scripts/probe_scooby_epicardioids_checkpoint_decoder.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/audit" \
  > "${STAGE}/validation/audit.stdout.json" \
  2> "${STAGE}/validation/audit.stderr.txt"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/probe_scooby_epicardioids_checkpoint_decoder.py" \
  "${ROOT}/tests/unit/test_probe_scooby_epicardioids_checkpoint_decoder.py" \
  "${ROOT}/slurm/probe_scooby_epicardioids_checkpoint_decoder_cpu.sbatch" \
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
    receipt["status"] != "pass_checkpoint_decoder_executable_full_sequence_and_native_context_still_blocked"
    or receipt["checkpoint"]["sha256"] != "f945e1937a46ddd5cb31249ef71e721e7a7b1381aa207fb2de449f23ac75a9e8"
    or not probe["strict_restore"]
    or probe["missing_keys"]
    or probe["unexpected_keys"]
    or probe["state_dict_tensors"] != 225
    or probe["state_dict_elements_including_buffers"] != 193241454
    or not probe["deterministic_repeat_bit_identical"]
    or not probe["checkpoint_deserialized_with_safetensors"]
    or not probe["checkpoint_decoder_forward_executed"]
    or probe["full_sequence_forward_executed"]
    or probe["native_cell_context_used"]
    or probe["project_data_read"]
    or probe["sealed_data_read"]
    or receipt["full_sequence_forward_executed"]
    or receipt["native_cell_context_used"]
    or receipt["project_data_read"]
    or receipt["sealed_data_read"]
    or receipt["task_scope"]["rna_conditioned_atac_eligible"]
    or receipt["license"]["open_champion_eligible"]
):
    raise SystemExit("Scooby checkpoint decoder receipt differs")
freeze_tree(
    stage,
    {
        "artifact_class": "restricted_observed_atac_model_checkpoint_decoder_probe",
        "model_id": "scooby_epicardioids",
        "checkpoint_sha256": receipt["checkpoint"]["sha256"],
        "checkpoint_deserialized_with_safetensors": True,
        "checkpoint_decoder_forward_executed": True,
        "full_sequence_forward_executed": False,
        "native_cell_context_used": False,
        "project_data_read": False,
        "sealed_data_read": False,
        "rna_conditioned_atac_eligible": False,
        "open_champion_eligible": False,
        "status": "pass_checkpoint_decoder_executable_full_sequence_and_native_context_still_blocked",
    },
)
publish_directory_noreplace(stage, target)
PY

trap - EXIT
printf 'Scooby strict checkpoint decoder probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
