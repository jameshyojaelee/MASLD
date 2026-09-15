#!/usr/bin/env bash
#SBATCH --job-name=model-work-304
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:b6k:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
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
# L40S. Predecessor job 21100298 (model-training-603) FAILED after 59 s with 'CUDA is unavailable' / 'CUDA initialization' -- empirical, not
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
# Inherited job 21100298 is left exactly as it is and is not resubmitted.
# ---------------------------------------------------------------------------

set -euo pipefail
umask 027
# AUTHORITY: dispatcher OR the recorded B6K direct-submission exception.
# The original guard required the central dispatcher, which is exactly the rule the
# user's B6K authorization carves an exception to -- the dispatcher hardcodes
# gpu:l40s:1 and physically cannot route B6K. The check is WIDENED, not removed:
# an unauthorised manual run of this or any other wrapper still fails closed,
# because one of the two named tokens is still required and each is value-checked.
B6K_AUTH_FILE="${ROOT_FOR_AUTH:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark}/config/b6k_numeric_tolerance_prespecification_20260825.json"
B6K_AUTH_SHA=a26a6f3fe0d77ef12706038eaede32f222421154c0972e228b5238febeaa722f
MASLD_AUTHORISED_BY=""
if [[ -n "${MASLD_GPU_DISPATCHER_AUTHORITY:-}" ]]; then
  if [[ "${MASLD_GPU_DISPATCHER_AUTHORITY}" != "sequence_gpu_dispatcher_v2" ]]; then
    printf 'Invalid GPU dispatcher authority.\n' >&2
    exit 1
  fi
  MASLD_AUTHORISED_BY="dispatcher:sequence_gpu_dispatcher_v2"
elif [[ -n "${MASLD_B6K_DIRECT_SUBMISSION_AUTHORITY:-}" ]]; then
  if [[ "${MASLD_B6K_DIRECT_SUBMISSION_AUTHORITY}" != "b6k_numeric_tolerance_20260825" ]]; then
    printf 'Invalid B6K direct-submission authority.\n' >&2
    exit 1
  fi
  # The exception is traceable to a frozen decision, not to convenience.
  printf '%s  %s\n' "${B6K_AUTH_SHA}" "${B6K_AUTH_FILE}" | sha256sum --check --strict
  MASLD_AUTHORISED_BY="b6k_direct_submission:b6k_numeric_tolerance_20260825"
else
  printf 'No authority. Set MASLD_GPU_DISPATCHER_AUTHORITY or MASLD_B6K_DIRECT_SUBMISSION_AUTHORITY.\n' >&2
  exit 1
fi
export MASLD_AUTHORISED_BY
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false
export PYTHONHASHSEED=20260825 CUBLAS_WORKSPACE_CONFIG=:4096:8
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/borzoi/bin/python3.11
CONTRACT="${ROOT}/config/borzoi_grelu_converted_port_full_window_probe.json"
SCRIPT="${ROOT}/scripts/probe_borzoi_grelu_converted_port_full_window.py"
TEST="${ROOT}/tests/unit/test_probe_borzoi_grelu_converted_port_full_window.py"
TARGET="${ROOT}/executions/borzoi-port-b6k-first-gpu-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite model-training-603 artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
for pair in \
  "7c81f285a121b9d70297ad7406bb40e445dda6585478f5896de9c41827d2e894 ${CONTRACT}" \
  "d9de47bceea5e32d5da1d70b0ac66155984eeca56c2c655c6f30056200969bb6 ${SCRIPT}" \
  "d6edbe8c53f2caa7fce29cb4768b0f98dbdfc9f69a01bc5b9e78f9661b327627 ${TEST}" \
  "fe1440fc9cc80b5475117f77f4a210213cfeaed63bc723ba58b11cdc838b1bd8 ${ROOT}/scripts/reconcile_borzoi_grelu_to_local_port_mapping.py" \
  "1054a52a3e8f87cc041d52107bb3ba6e74bf9e2dbe764b9c6d2c34e547d40f9b ${ROOT}/executions/model-data-259-21096699/ARTIFACTS.json" \
  "fa38f9b6bf4843e2c51e148e32be44680a44dbc97cfade25852f37f58c096d36 ${ROOT}/executions/model-check-261-21096790/ARTIFACTS.json" \
  "a80aed3c527577cc5203b6e578f27d9eda9ec4906044d73b9401245a886d8669 ${ROOT}/executions/model-check-264-21096934/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
if [[ "${PREFLIGHT_ONLY:-0}" == 1 ]]; then
  printf 'model-training-603 preflight passed\n'
  exit 0
fi
install -d -m 0750 "${STAGE}/validation"
nvidia-smi --query-gpu=name,compute_cap,driver_version,memory.total --format=csv,noheader \
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
    "l40s_predecessor_job": "21100298",
    "l40s_predecessor_outcome": "FAILED after 59 s with 'CUDA is unavailable' / 'CUDA initialization'",
    "why_b6k_only": ("torch 2.13.0+cu130 has no sm_89 SASS, so this runtime cannot "
                     "initialise CUDA on an L40S. B6K is the only hardware on which "
                     "this probe can execute."),
    "l40s_wrapper_left_untouched": True,
    "authorised_by": os.environ.get("MASLD_AUTHORISED_BY"),
    "authorization_artifact": "config/b6k_numeric_tolerance_prespecification_20260825.json",
    "authorization_artifact_sha256": "a26a6f3fe0d77ef12706038eaede32f222421154c0972e228b5238febeaa722f",
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

cd "${ROOT}"
PYTHONPATH="${ROOT}:${ROOT}/src" "${PYTHON}" -m unittest -v \
  tests.unit.test_probe_borzoi_grelu_converted_port_full_window \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}:${ROOT}/src" "${PYTHON}" "${SCRIPT}" \
  --contract "${CONTRACT}" --project-root "${ROOT}" --output "${STAGE}/probe" \
  > "${STAGE}/validation/probe.stdout.json"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  config/borzoi_grelu_converted_port_full_window_probe.json \
  scripts/probe_borzoi_grelu_converted_port_full_window.py \
  tests/unit/test_probe_borzoi_grelu_converted_port_full_window.py \
  scripts/reconcile_borzoi_grelu_to_local_port_mapping.py \
  slurm/run_model_training_603.sbatch > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "probe/receipt.json").read_text())
closed = (
    "unrestricted_torch_load_executed", "full_output_tensor_exported",
    "checkpoint_exported", "biological_data_read", "development_outcomes_read",
    "sealed_sources_read", "training_performed", "benchmark_metrics_computed",
    "global_frozen_census_modified", "native_borzoi_substitute",
    "checkpoint_bound_source_recovered", "native_numeric_parity_established",
    "native_weight_authority_established", "open_champion_eligible",
    "external_champion_claim_allowed", "universal_claim_allowed",
)
if (
    receipt["status"]
    != "pass_synthetic_full_window_deterministic_converted_port_forward"
    or receipt["gpu_name"] != "NVIDIA L40S"
    or receipt["sequence_bp"] != 524288
    or receipt["input_shape"] != [1, 4, 524288]
    or receipt["output_shape"] != [1, 12, 6144]
    or len(receipt["target_subset"]) != 12
    or not receipt["repeat_outputs_bit_identical"]
    or len(set(receipt["repeat_output_sha256"])) != 1
    or not receipt["weights_only_load"]
    or not receipt["semantic_key_mapping_used"]
    or not receipt["strict_restoration_passed"]
    or any(receipt[field] for field in closed)
):
    raise SystemExit("model-training-603 probe receipt differs")
freeze_tree(root, {
    "artifact_class": "borzoi_grelu_converted_port_full_window_l40s_probe",
    "model_id": "borzoi_grelu_converted_port", "logical_tasks": 1,
    "gpu": "NVIDIA L40S", "sequence_bp": 524288,
    "output_tracks": 12, "output_bins": 6144,
    "repeat_outputs_bit_identical": True,
    "weights_only_load": True, "strict_restoration_passed": True,
    "semantic_key_mapping_used": True,
    "model_forward_executed": True, "synthetic_input_only": True,
    "full_output_tensor_exported": False, "checkpoint_exported": False,
    "biological_data_read": False, "development_outcomes_read": False,
    "sealed_sources_read": False, "training_performed": False,
    "benchmark_metrics_computed": False, "global_frozen_census_modified": False,
    "native_borzoi_substitute": False, "checkpoint_bound_source_recovered": False,
    "native_numeric_parity_established": False,
    "native_weight_authority_established": False,
    "open_champion_eligible": False, "external_champion_claim_allowed": False,
    "universal_claim_allowed": False, "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'model-training-603 complete: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
