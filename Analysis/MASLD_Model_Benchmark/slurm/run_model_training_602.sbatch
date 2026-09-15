#!/usr/bin/env bash
#SBATCH --job-name=model-training-602
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
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/corgi-nvidia-pytorch-24.01-py3-recovery-21064690"
IMAGE="${RUNTIME}/nvidia-pytorch-24.01-py3.sif"
OVERLAY="${ROOT}/executions/environments/corgi-rotary-emb-0.1-flashattn2.0.4-build-r1-21064845"
ACQUISITION="${ROOT}/executions/upstream-acquisitions/corgi-d84c000-zenodo18630048-21064463"
SOURCE="${ACQUISITION}/source/corgi-d84c00082da3d04eb05169e290e53f1c7da2b9c5"
CHECKPOINT="${ACQUISITION}/downloads/corgi_model.pt"
FIXTURE="${ROOT}/executions/data-prep-264-21096833"
TILE_CONTRACT="${ROOT}/executions/corgi-outcome-aligned-tiles-21068311"
TILE_ROSTER="${ROOT}/executions/corgi-tile-smoke-subset-21068454"
ADAPTER="${ROOT}/src/masld_bench/adapters/corgi_film_head.py"
TARGET="${ROOT}/executions/model-training-602-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite model-training-602 artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
for pair in \
  "01de778c17c67d92de9a2e242350e17443dbfe76b7c5a5f8bc447eae5408cd3c ${RUNTIME}/ARTIFACTS.json" \
  "eace72d1b831a4834d692134749bb82093c29930fb3a0258ae667b8acb840959 ${IMAGE}" \
  "c60c4ebf0b32df7d158a3a1d781f5d45fdf4187ab0fba9f4bea35605a92e4c3a ${OVERLAY}/ARTIFACTS.json" \
  "73aad9ab9ad3e080a2ccaa0cee1a1d727d126a454f0e69f74f343a900b168c0f ${ACQUISITION}/ARTIFACTS.json" \
  "09d8811b612cba4f75361feeba29bd0676efde83f0e4785a37990c2d2f52a0a2 ${ROOT}/executions/corgi-native-runtime-l40s-r8-21065022/ARTIFACTS.json" \
  "ee695f9e421c140d77caf807de8fc05739bff4bbeac8a045b4cbdc80fa21ae7c ${ROOT}/executions/model-check-262-21096805/ARTIFACTS.json" \
  "63700ed1ad974175c65322047bde7ac39f067971ca25f90d33427ea0b5a1d20a ${FIXTURE}/ARTIFACTS.json" \
  "1cdf94aa4c37ead139125ffc3d8fee358605a2394a206b1f3e79190ecd25645c ${TILE_CONTRACT}/ARTIFACTS.json" \
  "2617f3a85a9ced8b0f7cb9efe3819b1e57177629221a2013e39fc7c96e187f5b ${TILE_ROSTER}/ARTIFACTS.json" \
  "6dd6242696f458efbec72bf73a36938802487a6d0e8f4e1ae2a401a3872129ca ${ADAPTER}" \
  "951f1349043dd9f8c2c5697bd7f64d93fa6f16b1f28b6dc337b1ef420d82f904 ${ROOT}/scripts/probe_corgi_regular_film_head_l40s.py" \
  "2259363a1f8d0729f095c904f1f59ced21c349aa8ec140b98a1445f30cf91dfd ${ROOT}/tests/unit/test_probe_corgi_regular_film_head_l40s.py" \
  "368585b968f690646ccd60a59df146c24aaf429b28058225491c1c991d3fe596 ${ROOT}/tests/unit/test_corgi_film_head_adapter.py" \
  "2a68bc359f246f4fecc6e835fc6fd798861b4c631fd9d665fccd85d200b6e31c ${ROOT}/scripts/corgi_predict_outcome_aligned_tiles.py" \
  "e351024065408c644b0476fffa41a002e6fd57f7b768972d020b3275fcfcf9f5 ${ROOT}/scripts/corgi_predict_valid_contexts.py" \
  "f21d31c0c9888e3a80f730ff57ba4ccaf6b57cad887d542d16a442a241f34be6 ${ROOT}/scripts/corgi_probe_native_runtime.py"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
if [[ "${PREFLIGHT_ONLY:-0}" == 1 ]]; then
  printf 'model-training-602 preflight passed\n'
  exit 0
fi
install -d -m 0750 "${STAGE}/validation"
module purge
module load singularity/4.2.2
nvidia-smi --query-gpu=name,compute_cap,driver_version,memory.total --format=csv,noheader \
  > "${STAGE}/validation/nvidia-smi.csv"
COMMON=(--nv --cleanenv --containall --net --network none --bind "${ROOT}:${ROOT}:ro" --bind "${STAGE}:${STAGE}:rw" --pwd "${ROOT}" --env "PYTHONPATH=${ROOT}:${OVERLAY}/site-packages" --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 --env PYTHONHASHSEED=20260825 --env CUBLAS_WORKSPACE_CONFIG=:4096:8)
singularity exec "${COMMON[@]}" "${IMAGE}" python -m unittest -v \
  tests.unit.test_corgi_film_head_adapter \
  tests.unit.test_probe_corgi_regular_film_head_l40s \
  > "${STAGE}/validation/tests.stdout.txt" 2> "${STAGE}/validation/tests.stderr.txt"
singularity exec "${COMMON[@]}" "${IMAGE}" python \
  scripts/probe_corgi_regular_film_head_l40s.py \
  --source "${SOURCE}" --checkpoint "${CHECKPOINT}" --adapter "${ADAPTER}" \
  --tile-contract "${TILE_CONTRACT}" --tile-roster "${TILE_ROSTER}" \
  --fixture "${FIXTURE}" --output "${STAGE}/probe" \
  > "${STAGE}/validation/probe.stdout.json"
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:ro" --pwd "${ROOT}" "${IMAGE}" \
  python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
cd "${ROOT}"
sha256sum src/masld_bench/adapters/corgi_film_head.py \
  scripts/probe_corgi_regular_film_head_l40s.py \
  tests/unit/test_corgi_film_head_adapter.py \
  tests/unit/test_probe_corgi_regular_film_head_l40s.py \
  slurm/run_model_training_602.sbatch > "${STAGE}/source.sha256"
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
    "development_atac_outcomes_read", "histology_or_disease_labels_read",
    "test_or_sealed_features_or_outcomes_read", "benchmark_metrics_computed",
    "production_training_performed", "production_prediction_performed",
    "promotion_gate_evaluated", "native_numeric_parity_established",
    "open_champion_eligible", "global_census_modified",
    "global_promotion_gate_modified",
)
if (
    receipt["status"] != "pass_full_window_single_step_and_resume"
    or receipt["gpu_name"] != "NVIDIA L40S"
    or receipt["sequence_bp"] != 524288
    or receipt["output_shape"] != [1, 1, 6144]
    or receipt["trainable_parameter_numel"] != 5869229
    or not receipt["batchnorm_running_statistics_frozen"]
    or receipt["frozen_parameter_gradients_present"]
    or not receipt["representative_frozen_parameters_unchanged"]
    or receipt["synthetic_probe_delta_reusable_as_benchmark_model"]
    or any(receipt[field] for field in closed)
):
    raise SystemExit("model-training-602 probe receipt differs")
freeze_tree(root, {
    "artifact_class": "corgi_regular_film_plus_head_l40s_full_window_probe",
    "model_id": "corgi_regular", "task_id": "rna_conditioned_atac",
    "adaptation_rung": "film_plus_head", "gpu": "NVIDIA L40S",
    "logical_tasks": 1, "sequence_bp": 524288,
    "trainable_parameter_numel": receipt["trainable_parameter_numel"],
    "released_row_initialization_parity_passed": True,
    "single_training_step_passed": True, "checkpoint_resume_passed": True,
    "synthetic_probe_delta_reusable_as_benchmark_model": False,
    "development_atac_outcomes_read": False,
    "histology_or_disease_labels_read": False,
    "test_or_sealed_features_or_outcomes_read": False,
    "benchmark_metrics_computed": False, "production_training_performed": False,
    "production_prediction_performed": False, "promotion_gate_evaluated": False,
    "native_numeric_parity_established": False, "open_champion_eligible": False,
    "global_census_modified": False, "global_promotion_gate_modified": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'model-training-602 complete: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
