#!/usr/bin/env bash
#SBATCH --job-name=masld-peakvi-context-smoke
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PREPARED="${ROOT}/executions/multivi-smoke-prepared-21065093"
PREPARED_ARTIFACTS_SHA256=6c10206033545fe7ed1ecbd9f7d64e107e2b3742de06db55a838bbdb46a1cdaa
RUNTIME="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
RUNTIME_ARTIFACTS_SHA256=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
SITE_PACKAGES="${RUNTIME}/site-packages"
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
TORCH_ARTIFACTS_SHA256=d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3
TORCH_SITE_PACKAGES="${TORCH}/site-packages"
TARGET="${ROOT}/executions/peakvi-training-context-smoke-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed PeakVI context transfer: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite PeakVI context transfer: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for record in \
  "${PREPARED_ARTIFACTS_SHA256} ${PREPARED}/ARTIFACTS.json" \
  "${RUNTIME_ARTIFACTS_SHA256} ${RUNTIME}/ARTIFACTS.json" \
  "${TORCH_ARTIFACTS_SHA256} ${TORCH}/ARTIFACTS.json"; do
  printf '%s\n' "${record}" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${PREPARED}" "${RUNTIME}" "${TORCH}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for path in sys.argv[1:]:
    verify_frozen_tree(Path(path))
PY
module purge
module load python/3.12.3-GCCcore-13.3.0
PYTHONPATH="${TORCH_SITE_PACKAGES}:${SITE_PACKAGES}:${ROOT}" python -m unittest -v \
  tests.unit.test_peakvi_context_transfer_smoke \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${TORCH_SITE_PACKAGES}:${SITE_PACKAGES}:${ROOT}/scripts" python \
  "${ROOT}/scripts/peakvi_context_transfer_smoke.py" \
  --prepared "${PREPARED}/prepared" \
  --output "${STAGE}/predictions" \
  --seed 1103 \
  --epochs 50 \
  > "${STAGE}/validation/fit_predict.stdout.json"
PYTHONPATH="${TORCH_SITE_PACKAGES}:${SITE_PACKAGES}" python -m pip freeze --all \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/peakvi_context_transfer_smoke.py" \
  "${ROOT}/scripts/multivi_fit_predict_smoke.py" \
  "${ROOT}/tests/unit/test_peakvi_context_transfer_smoke.py" \
  "${ROOT}/slurm/peakvi_context_transfer_smoke.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench import PredictionBundle
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "predictions/receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["model_id"] != "peakvi_training_context"
    or receipt["outer_unit"] != "donor"
    or set(receipt["folds"]) != set(map(str, range(5)))
    or receipt["query_rna_used"]
    or receipt["held_atac_available_to_adapter"]
    or receipt["outcomes_read"]
    or receipt["metrics_calculated_by_adapter"]
    or receipt["rna_conditioned_atac_eligible"]
    or receipt["sealed_inference_eligible"]
    or receipt["champion_claim_allowed"]
):
    raise SystemExit("PeakVI context-transfer receipt differs")
for fold, record in receipt["folds"].items():
    if record["held_atac_read"] or not record["repeat_prototype_decode_bit_identical"]:
        raise SystemExit(f"PeakVI fold receipt differs: {fold}")
    PredictionBundle.load_json(root / "predictions" / record["prediction_bundle"])
freeze_tree(root, {
    "artifact_class": "peakvi_training_context_smoke",
    "model_id": "peakvi_training_context",
    "dataset_id": "gse296875",
    "prepared_artifacts_sha256": "6c10206033545fe7ed1ecbd9f7d64e107e2b3742de06db55a838bbdb46a1cdaa",
    "held_atac_available_to_adapter": False,
    "outcomes_read": False,
    "rna_conditioned_atac_eligible": False,
    "sealed_inference_eligible": False,
    "champion_claim_allowed": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'PeakVI context-transfer artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
