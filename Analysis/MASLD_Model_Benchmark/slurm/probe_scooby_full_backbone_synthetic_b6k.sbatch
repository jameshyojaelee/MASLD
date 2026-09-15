#!/usr/bin/env bash
#SBATCH --job-name=model-probe-016
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:b6k:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err


# ---------------------------------------------------------------------------
# B6K MIGRATION. Separate wrapper by necessity, not style: the CPU queue
# validators assert that their L40S target contains the literal
# the L40S gres directive verbatim, so an in-place edit would break the target
# (the literal is deliberately not reproduced here: a wrapper quoting it would
#  itself register as L40S-targeting under a naive grep)
# wrapper's own validator. The L40S wrapper is left untouched.
#
# WHY B6K IS THE ONLY HARDWARE THIS CAN RUN ON:
# the interpreter is micromamba envs/borzoi at torch 2.13.0+cu130, whose arch
# list is sm_75/80/86/90/100/120 with NO sm_89. It cannot initialise CUDA on an
# L40S. Predecessor job 21100299 (model-probe-284) FAILED after 10 s on ne1dg6-009 -- empirical, not
# inferred. This is one of exactly three wrappers in the repository using a
# cu130 runtime bare; every other cu130 consumer layers the
# pytorch-2.6.0-cu124-portable overlay and is L40S-safe.
#
# FIRST GPU EXECUTION of this runtime on B6K. Behaviour is being established,
# not re-established across a change.
#
# Carries no Resource separation baseline (zero references to resource_snapshot /
# assert_resource_unchanged / candidates/), so the active Fig4h Resource
# adoption cannot invalidate it.
#
# Direct sbatch; never enters the dispatcher roster, which hardcodes gpu:l40s:1.
# Inherited job 21100299 is left exactly as it is and is not resubmitted.
# ---------------------------------------------------------------------------

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
# B6K requirements variant: the L40S original is unmodified and stays hash-pinned by
# scripts/audit_scooby_biological_fixture_readiness.py. required_gpu_name_contains
# remains a REQUIRED assertion, retargeted rather than relaxed.
CONFIG="${ROOT}/config/scooby_epicardioids_full_backbone_synthetic_probe_b6k.json"
TARGET="${ROOT}/executions/scooby-b6k-first-gpu-${SLURM_JOB_ID}"
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

# Layer guard: prove torch resolves inside the intended runtime and that a real
# Blackwell device is live, before any scientific work runs.
"${PYTHON}" - "${STAGE}" <<'LAYERGUARD'
import json, os, sys
from pathlib import Path
import torch
stage = Path(sys.argv[1])
resolved = os.path.realpath(torch.__file__)
g = {
    "arch_list": list(torch.cuda.get_arch_list()),
    "cuda_available": torch.cuda.is_available(),
    "cuda_runtime": torch.version.cuda,
    "expected_runtime_substring": "envs/borzoi",
    "layer_guard_passed": "envs/borzoi" in resolved,
    "torch_file": resolved,
    "torch_version": torch.__version__,
}
if torch.cuda.is_available():
    g["device_name"] = torch.cuda.get_device_name(0)
    g["device_capability"] = list(torch.cuda.get_device_capability(0))
    g["device_count"] = torch.cuda.device_count()
g["b6k_migration"] = {
    "first_gpu_execution_of_this_runtime": True,
    "l40s_predecessor_job": "21100299",
    "l40s_predecessor_outcome": "FAILED after 10 s on ne1dg6-009",
    "why_b6k_only": ("torch 2.13.0+cu130 has no sm_89 SASS, so this runtime cannot "
                     "initialise CUDA on an L40S. B6K is the only hardware on which "
                     "this probe can execute."),
    "l40s_wrapper_left_untouched": True,
}
(stage / "validation").mkdir(parents=True, exist_ok=True)
(stage / "validation" / "b6k_layer_guard.json").write_text(
    json.dumps(g, indent=2, sort_keys=True) + "\n")
print(json.dumps(g, indent=2, sort_keys=True))
if not g["layer_guard_passed"]:
    raise SystemExit(f"layer guard failed: torch resolved to {resolved}")
if torch.__version__ != "2.13.0+cu130":
    raise SystemExit(f"unexpected torch: {torch.__version__}")
if not g["cuda_available"] or g.get("device_count") != 1:
    raise SystemExit("exactly one visible CUDA device is required")
if g["device_capability"][0] != 12:
    raise SystemExit(f"expected Blackwell sm_12x, got {g.get('device_name')}")
if "sm_89" in g["arch_list"]:
    raise SystemExit("unexpected sm_89 SASS; this build was recorded as L40S-incapable")
LAYERGUARD

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
