#!/usr/bin/env bash
#SBATCH --job-name=model-probe-284
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
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
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/borzoi/bin/python
CONFIG="${ROOT}/config/scooby_epicardioids_full_backbone_synthetic_probe.json"
TARGET="${ROOT}/executions/model-probe-284-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite full-backbone probe: %s\n' "${TARGET}" >&2
  exit 1
fi
for pair in \
  "657f20fb7fe110fdd4f690b874d450f73fc6fc9b3e40c0276f0bade73bc9047a ${ROOT}/executions/model-data-281-21097353/ARTIFACTS.json" \
  "17aad708f568157785905fbed57dc78fc6bf0af01f7aacaa493108340615a88a ${ROOT}/executions/model-check-282-21097437/ARTIFACTS.json" \
  "835244d33a0fcc2808efc76551eb0f0f771e08ec17e0acf12599cd7ceb9365ea ${ROOT}/executions/model-check-279-21097163/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" - "${ROOT}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree


root = Path(sys.argv[1])
for value in (
    "executions/model-data-281-21097353",
    "executions/model-check-282-21097437",
    "executions/model-check-279-21097163",
):
    verify_frozen_tree(root / value)
PY
if [[ "${PREFLIGHT_ONLY:-0}" == 1 ]]; then
  printf 'Scooby full-backbone preflight passed\n'
  exit 0
fi

install -d -m 0750 "${STAGE}/validation"
nvidia-smi --query-gpu=name,uuid,driver_version,memory.total,compute_cap --format=csv,noheader \
  > "${STAGE}/validation/nvidia-smi.csv"
"${PYTHON}" -m unittest -v \
  tests.unit.test_probe_scooby_full_backbone_synthetic \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${PYTHON}" "${ROOT}/scripts/probe_scooby_full_backbone_synthetic.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/audit" \
  > "${STAGE}/validation/audit.stdout.json" \
  2> "${STAGE}/validation/audit.stderr.txt"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/probe_scooby_full_backbone_synthetic.py" \
  "${ROOT}/tests/unit/test_probe_scooby_full_backbone_synthetic.py" \
  "${ROOT}/slurm/probe_scooby_full_backbone_synthetic_l40s.sbatch" \
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
    receipt["status"] != "pass_full_backbone_synthetic_native_liver_context_still_blocked"
    or probe["output"]["shape"] != [1, 6144, 3]
    or not probe["output"]["finite"]
    or not probe["deterministic_repeat_bit_identical"]
    or not probe["strict_restore"]
    or probe["missing_keys"]
    or probe["unexpected_keys"]
    or not probe["checkpoint_deserialized_with_safetensors"]
    or not probe["full_sequence_forward_executed"]
    or probe["native_cell_context_used"]
    or probe["biological_task_executed"]
    or probe["model_fitted_or_adapted"]
    or probe["project_data_read"]
    or probe["sealed_data_read"]
    or receipt["native_cell_context_used"]
    or receipt["biological_task_executed"]
    or receipt["model_fitted_or_adapted"]
    or receipt["project_data_read"]
    or receipt["sealed_data_read"]
    or receipt["task_scope"]["rna_conditioned_atac_eligible"]
):
    raise SystemExit("Scooby full-backbone receipt differs")
freeze_tree(
    stage,
    {
        "artifact_class": "restricted_observed_atac_model_full_backbone_synthetic_probe",
        "model_id": "scooby_epicardioids",
        "checkpoint_sha256": receipt["checkpoint"]["sha256"],
        "output_sha256": probe["output"]["sha256"],
        "full_sequence_forward_executed": True,
        "deterministic_repeat_bit_identical": True,
        "native_cell_context_used": False,
        "biological_task_executed": False,
        "model_fitted_or_adapted": False,
        "project_data_read": False,
        "sealed_data_read": False,
        "rna_conditioned_atac_eligible": False,
        "open_champion_eligible": False,
        "status": "pass_full_backbone_synthetic_native_liver_context_still_blocked",
    },
)
publish_directory_noreplace(stage, target)
PY

trap - EXIT
printf 'Scooby full-backbone synthetic probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
