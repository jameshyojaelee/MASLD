#!/usr/bin/env bash
#SBATCH --job-name=masld-midas-fold
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2
export CUBLAS_WORKSPACE_CONFIG=:4096:8

: "${MIDAS_FOLD:?MIDAS_FOLD is required}"
: "${MIDAS_SEED:?MIDAS_SEED is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PREPARED="${ROOT}/executions/multivi-smoke-prepared-21065093"
PREPARED_ARTIFACTS_SHA256=6c10206033545fe7ed1ecbd9f7d64e107e2b3742de06db55a838bbdb46a1cdaa
MIDAS="${ROOT}/executions/environments/midas-v0.3.0-portable-21066061"
MIDAS_ARTIFACTS_SHA256=c861185715ed6952f8a3a09b37b04e17e5a2e8e22a23a2e3c8ab0688508fddce
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
BASE_ARTIFACTS_SHA256=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
TORCH_ARTIFACTS_SHA256=d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3
SITE_PACKAGES="${MIDAS}/site-packages:${TORCH}/site-packages:${BASE}/site-packages"
TARGET="${ROOT}/executions/midas-inductive-seed${MIDAS_SEED}-fold${MIDAS_FOLD}-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed MIDAS donor fold: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite MIDAS donor fold: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for record in \
  "${PREPARED_ARTIFACTS_SHA256} ${PREPARED}/ARTIFACTS.json" \
  "${MIDAS_ARTIFACTS_SHA256} ${MIDAS}/ARTIFACTS.json" \
  "${BASE_ARTIFACTS_SHA256} ${BASE}/ARTIFACTS.json" \
  "${TORCH_ARTIFACTS_SHA256} ${TORCH}/ARTIFACTS.json"; do
  printf '%s\n' "${record}" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${PREPARED}" "${MIDAS}" "${BASE}" "${TORCH}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for path in sys.argv[1:]:
    verify_frozen_tree(Path(path))
PY
module purge
module load python/3.12.3-GCCcore-13.3.0
PYTHONPATH="${SITE_PACKAGES}:${ROOT}" python -m unittest -v \
  tests.unit.test_midas_fit_predict_fold \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${SITE_PACKAGES}:${ROOT}/scripts" python \
  "${ROOT}/scripts/midas_fit_predict_fold.py" \
  --prepared "${PREPARED}/prepared" \
  --helper-path "${ROOT}/scripts/multivi_fit_predict_smoke.py" \
  --output "${STAGE}/predictions" \
  --fold "${MIDAS_FOLD}" --seed "${MIDAS_SEED}" --epochs 30 \
  > "${STAGE}/validation/fit_predict.stdout.json" \
  2> "${STAGE}/validation/fit_predict.stderr.txt"
PYTHONPATH="${SITE_PACKAGES}" python -m pip freeze --all \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/midas_fit_predict_fold.py" \
  "${ROOT}/scripts/multivi_fit_predict_smoke.py" \
  "${ROOT}/tests/unit/test_midas_fit_predict_fold.py" \
  "${ROOT}/slurm/midas_fit_predict_fold.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${MIDAS_FOLD}" "${MIDAS_SEED}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench import PredictionBundle
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
fold = int(sys.argv[2])
seed = int(sys.argv[3])
receipt = json.loads((root / "predictions/receipt.json").read_text())
record = receipt["folds"].get(str(fold), {})
if (
    receipt["status"] != "pass"
    or receipt["model_id"] != "midas_inductive"
    or receipt["exact_distribution"] != "scmidas 0.3.0"
    or receipt["pairing_topology"] != "same_nucleus_training"
    or receipt["outer_unit"] != "donor"
    or set(receipt["folds"]) != {str(fold)}
    or receipt["base_seed"] != seed
    or not receipt["official_lightning_trainer_used"]
    or not receipt["training_datalist_paired_only"]
    or receipt["query_rna_used_for_training"]
    or receipt["query_atac_state"] != "structurally_missing"
    or receipt["held_atac_input_exposed"]
    or receipt["outcomes_read"]
    or receipt["metrics_calculated_by_adapter"]
    or receipt["test_outcomes_read"]
    or receipt["sealed_inference_eligible"]
    or receipt["champion_claim_allowed"]
    or not record.get("strict_state_resume")
    or not record.get("repeated_query_prediction_bit_identical")
    or record.get("query_rna_used_for_training")
    or record.get("held_atac_read")
):
    raise SystemExit("MIDAS donor-held fold receipt differs")
PredictionBundle.load_json(root / "predictions" / record["prediction_bundle"])
freeze_tree(root, {
    "artifact_class": "midas_inductive_fit_predict_fold",
    "model_id": "midas_inductive",
    "dataset_id": "gse296875",
    "fold": fold,
    "seed": seed,
    "prepared_artifacts_sha256": "6c10206033545fe7ed1ecbd9f7d64e107e2b3742de06db55a838bbdb46a1cdaa",
    "runtime_artifacts_sha256": "c861185715ed6952f8a3a09b37b04e17e5a2e8e22a23a2e3c8ab0688508fddce",
    "query_rna_used_for_training": False,
    "held_atac_input_exposed": False,
    "outcomes_read": False,
    "test_outcomes_read": False,
    "sealed_inference_eligible": False,
    "champion_claim_allowed": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'MIDAS donor-held fold: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
