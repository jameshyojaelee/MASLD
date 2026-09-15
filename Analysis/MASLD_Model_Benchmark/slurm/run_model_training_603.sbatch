#!/usr/bin/env bash
#SBATCH --job-name=model-training-603
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
: "${MASLD_GPU_DISPATCHER_AUTHORITY:?dispatcher authority is required}"
if [[ "${MASLD_GPU_DISPATCHER_AUTHORITY}" != "sequence_gpu_dispatcher_v2" ]]; then
  printf 'Invalid GPU dispatcher authority.\n' >&2
  exit 1
fi
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false
export PYTHONHASHSEED=20260825 CUBLAS_WORKSPACE_CONFIG=:4096:8
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/borzoi/bin/python3.11
CONTRACT="${ROOT}/config/borzoi_grelu_converted_port_full_window_probe.json"
SCRIPT="${ROOT}/scripts/probe_borzoi_grelu_converted_port_full_window.py"
TEST="${ROOT}/tests/unit/test_probe_borzoi_grelu_converted_port_full_window.py"
TARGET="${ROOT}/executions/model-training-603-${SLURM_JOB_ID}"
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
