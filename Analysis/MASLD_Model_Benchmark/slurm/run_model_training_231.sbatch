#!/usr/bin/env bash
#SBATCH --job-name=model-training-231
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
: "${MASLD_GPU_DISPATCHER_AUTHORITY:?dispatcher authority is required}"
if [[ "${MASLD_GPU_DISPATCHER_AUTHORITY}" != "sequence_gpu_dispatcher_v2" ]]; then
  printf 'Invalid GPU dispatcher authority.\n' >&2
  exit 2
fi

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
CONVERSION="${ROOT}/executions/environments/dna-lm-torch2.6-conversion-21064945"
CONFIG="${ROOT}/config/gse281364_task_native_sequence_control_production.toml"
TASK="${ROOT}/executions/model-check-227-21092007"
OUTCOMES="${ROOT}/executions/gse281364-replicate-outcomes-21066307"
PREFLIGHT="${ROOT}/executions/model-check-232-21092211"
TARGET="${ROOT}/executions/model-training-231-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed sequence-control production attempt: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite sequence-control production output.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/parts" "${STAGE}/validation"

for binding in \
  "d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3 executions/environments/pytorch-2.6.0-cu124-portable-21064747" \
  "5fed22c975200495ea53b38b09f5286cf13a9230a138cb72a0c758c85d5cb9d4 executions/environments/dna-lm-torch2.6-conversion-21064945" \
  "30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af executions/environments/scvi-tools-1.5.0.post1-portable-21064579" \
  "622ca08fa8df12f91349431c40e82bde693c56423c755dea6f72dd76c5478a5f executions/model-check-227-21092007" \
  "cb42bd4037127a4e4d6a406d6f7b00c706bc93b631cdf93083ce994809e791db executions/gse281364-replicate-outcomes-21066307" \
  "24a98562653c30fe44dd0256460a923fbaf7cfd479ca1f12707dda83e4e09ceb executions/model-check-232-21092211"; do
  read -r expected relative <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}/ARTIFACTS.json" \
    | sha256sum --check --strict
done
for binding in \
  "845820b7d7a510e6b959bb318914da161165fd71771f0a139b74c31da96bafd9 config/gse281364_task_native_sequence_control_production.toml" \
  "58b83a584fab05cded6f254817548125a3fb489abe230e783c7b1d88cbf7cb15 scripts/gse281364_task_native_sequence_controls.py" \
  "79add1142cdaa7a5b4b4849d33e26b78ea91b91b50a7f065dd2296cbe470ffa9 scripts/materialize_gse281364_sequence_control_training_targets.py" \
  "9b37c5050d36c394bcb3c916d4743e6e5ea7ffe3178ab211146f460e85e6eb07 scripts/train_predict_gse281364_task_native_sequence_controls.py" \
  "7c274121e1fb63d1db1c0bf53b2efa3ac19c4674728317345c3278f545a78535 scripts/aggregate_gse281364_task_native_sequence_control_predictions.py" \
  "93eb6889e025d28c9da6bd522f52b6c3c04570d154a2d1f21ee3c2bc3a6c6a49 tests/unit/test_gse281364_task_native_sequence_controls.py" \
  "f7e5efb483021a2f644c35ae54d6f758e6693979664f324d3eaf2c6da2e6e831 tests/unit/test_train_predict_gse281364_task_native_sequence_controls.py" \
  "ed4bd726db8210d11c034b806c54fa541afce62664ee0846912022036e414635 scripts/offline_guard/sitecustomize.py"; do
  read -r expected relative <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/scripts/offline_guard:${CONVERSION}/site-packages:${TORCH}/site-packages:${BASE}/site-packages:${ROOT}/src:${ROOT}"
export XDG_CACHE_HOME="${SLURM_TMPDIR:-/tmp}/sequence-control-cache-${SLURM_JOB_ID}"
export TORCHINDUCTOR_CACHE_DIR="${XDG_CACHE_HOME}/torchinductor"
export TRITON_CACHE_DIR="${XDG_CACHE_HOME}/triton"
install -d -m 0700 "${XDG_CACHE_HOME}" "${TORCHINDUCTOR_CACHE_DIR}" "${TRITON_CACHE_DIR}"

python - <<'PY' > "${STAGE}/validation/gpu_runtime.json"
import json
import torch
if not torch.cuda.is_available():
    raise SystemExit("CUDA unavailable")
capability = torch.cuda.get_device_capability(0)
if torch.cuda.get_device_name(0).lower().find("l40s") < 0:
    raise SystemExit("allocated GPU is not L40S")
x = torch.arange(16, device="cuda", dtype=torch.float32).reshape(4, 4)
y = x @ x.T
torch.cuda.synchronize()
print(json.dumps({
    "status": "pass",
    "torch_version": torch.__version__,
    "cuda_version": torch.version.cuda,
    "device_name": torch.cuda.get_device_name(0),
    "compute_capability": list(capability),
    "finite_kernel": bool(torch.isfinite(y).all().item()),
}, indent=2, sort_keys=True))
PY
nvidia-smi -q > "${STAGE}/validation/nvidia-smi-q.txt"
python -m unittest -v \
  tests.unit.test_gse281364_task_native_sequence_controls \
  tests.unit.test_train_predict_gse281364_task_native_sequence_controls \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

PRIVATE="$(mktemp -d -p "${SLURM_TMPDIR:-/tmp}" sequence-target-production.XXXXXX)"
for fold in fold-0 fold-1 fold-2 fold-3 fold-4; do
  python "${ROOT}/scripts/materialize_gse281364_sequence_control_training_targets.py" \
    --project-root "${ROOT}" \
    --config "${CONFIG}" \
    --outer-fold "${fold}" \
    --output "${PRIVATE}/${fold}" \
    > "${STAGE}/validation/${fold}.materialize.stdout.json" \
    2> "${STAGE}/validation/${fold}.materialize.stderr.txt"
  python "${ROOT}/scripts/train_predict_gse281364_task_native_sequence_controls.py" \
    --config "${CONFIG}" \
    --frozen-task "${TASK}" \
    --training-targets "${PRIVATE}/${fold}" \
    --outer-fold "${fold}" \
    --device cuda \
    --output "${STAGE}/parts/${fold}" \
    > "${STAGE}/validation/${fold}.training.stdout.jsonl" \
    2> "${STAGE}/validation/${fold}.training.stderr.txt"
  find "${PRIVATE}/${fold}" -type f -delete
  find "${PRIVATE}/${fold}" -depth -type d -empty -delete
  if [[ -e "${PRIVATE}/${fold}" ]]; then
    printf 'Private target view remains after fold training: %s\n' "${fold}" >&2
    exit 1
  fi
done
rmdir "${PRIVATE}"
if [[ -e "${PRIVATE}" ]]; then
  printf 'Private target root remains after production cleanup\n' >&2
  exit 1
fi

python "${ROOT}/scripts/aggregate_gse281364_task_native_sequence_control_predictions.py" \
  --parts "${STAGE}/parts" \
  --output "${STAGE}/campaign" \
  > "${STAGE}/validation/aggregate.stdout.json" \
  2> "${STAGE}/validation/aggregate.stderr.txt"
python - "${STAGE}" <<'PY' > "${STAGE}/validation/final_contract.json"
from pathlib import Path
import json
import sys

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "campaign/receipt.json").read_text())
checkpoints = list((root / "parts").rglob("*.safetensors"))
fit_receipts = list((root / "parts").glob("*/training_receipts/*/seed-*.json"))
raw_targets = list(root.rglob("training_targets.tsv"))
if (
    receipt["status"] != "pass_prediction_only_production_campaign"
    or receipt["fit_count"] != 50
    or receipt["checkpoint_count"] != 50
    or receipt["training_receipt_count"] != 50
    or receipt["prediction_rows"] != 20660
    or len(checkpoints) != 50
    or len(fit_receipts) != 50
    or raw_targets
    or receipt["held_test_outcomes_visible_to_training_adapter"]
    or receipt["raw_training_targets_retained"]
    or receipt["outcomes_copied_to_output"]
    or receipt["benchmark_metrics_calculated"]
    or receipt["evaluator_invoked"]
    or receipt["model_ranking_performed"]
    or receipt["champion_eligible"]
    or receipt["external_claim_eligible"]
):
    raise SystemExit("final prediction-only contract differs")
print(json.dumps({
    "status": "pass",
    "fit_count": 50,
    "checkpoint_count": len(checkpoints),
    "training_receipt_count": len(fit_receipts),
    "prediction_rows": receipt["prediction_rows"],
    "raw_training_target_files": len(raw_targets),
    "benchmark_metrics_calculated": False,
}, indent=2, sort_keys=True))
PY
python -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/gse281364_task_native_sequence_controls.py" \
  "${ROOT}/scripts/materialize_gse281364_sequence_control_training_targets.py" \
  "${ROOT}/scripts/train_predict_gse281364_task_native_sequence_controls.py" \
  "${ROOT}/scripts/aggregate_gse281364_task_native_sequence_control_predictions.py" \
  "${ROOT}/tests/unit/test_gse281364_task_native_sequence_controls.py" \
  "${ROOT}/tests/unit/test_train_predict_gse281364_task_native_sequence_controls.py" \
  "${ROOT}/slurm/run_model_training_231.sbatch" \
  > "${STAGE}/source.sha256"
module purge

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "campaign/receipt.json").read_text())
freeze_tree(root, {
    "artifact_class": "gse281364_task_native_sequence_control_predictions",
    "dataset_id": "gse281364",
    "models": ["sequence_cnn_control", "sequence_transformer_control"],
    "fit_count": 50,
    "checkpoint_count": 50,
    "training_receipt_count": 50,
    "prediction_rows": 20660,
    "predictions_sha256": receipt["predictions_sha256"],
    "held_test_outcomes_visible_to_training_adapter": False,
    "raw_training_targets_retained": False,
    "benchmark_metrics_calculated": False,
    "evaluator_invoked": False,
    "champion_eligible": False,
    "external_claim_eligible": False,
    "status": "passed",
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Sequence-control production campaign: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
