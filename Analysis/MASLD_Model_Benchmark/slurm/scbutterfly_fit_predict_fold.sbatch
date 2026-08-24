#!/usr/bin/env bash
#SBATCH --job-name=masld-scbutterfly-fold
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export CUBLAS_WORKSPACE_CONFIG=:4096:8

: "${SCBUTTERFLY_FOLD:?SCBUTTERFLY_FOLD is required}"
: "${SCBUTTERFLY_SEED:?SCBUTTERFLY_SEED is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PREPARED="${ROOT}/executions/multivi-smoke-prepared-21065093"
PREPARED_ARTIFACTS_SHA256=6c10206033545fe7ed1ecbd9f7d64e107e2b3742de06db55a838bbdb46a1cdaa
SCBUTTERFLY="${ROOT}/executions/environments/scbutterfly-v0.0.9-portable-21066062"
SCBUTTERFLY_ARTIFACTS_SHA256=b495d6a28bff02387d4c1c56993df29bbac1cdc28119e8ec9bc4cff621965014
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
BASE_ARTIFACTS_SHA256=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
TORCH_ARTIFACTS_SHA256=d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3
SITE_PACKAGES="${SCBUTTERFLY}/site-packages:${TORCH}/site-packages:${BASE}/site-packages"
COMPONENT="${SCBUTTERFLY}/site-packages/scButterfly/model_component.py"
TARGET="${ROOT}/executions/scbutterfly-b-seed${SCBUTTERFLY_SEED}-fold${SCBUTTERFLY_FOLD}-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed scButterfly fold: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite scButterfly fold: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for record in \
  "${PREPARED_ARTIFACTS_SHA256} ${PREPARED}/ARTIFACTS.json" \
  "${SCBUTTERFLY_ARTIFACTS_SHA256} ${SCBUTTERFLY}/ARTIFACTS.json" \
  "${BASE_ARTIFACTS_SHA256} ${BASE}/ARTIFACTS.json" \
  "${TORCH_ARTIFACTS_SHA256} ${TORCH}/ARTIFACTS.json"; do
  printf '%s\n' "${record}" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${PREPARED}" "${SCBUTTERFLY}" "${BASE}" "${TORCH}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for path in sys.argv[1:]:
    verify_frozen_tree(Path(path))
PY
module purge
module load python/3.12.3-GCCcore-13.3.0
PYTHONPATH="${SITE_PACKAGES}:${ROOT}" python -m unittest -v \
  tests.unit.test_scbutterfly_fit_predict_smoke \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${SITE_PACKAGES}:${ROOT}/scripts" python \
  "${ROOT}/scripts/scbutterfly_fit_predict_smoke.py" \
  --component-path "${COMPONENT}" --self-test \
  > "${STAGE}/validation/self_test.stdout.txt" \
  2> "${STAGE}/validation/self_test.stderr.txt"
PYTHONPATH="${SITE_PACKAGES}:${ROOT}/scripts" python \
  "${ROOT}/scripts/scbutterfly_fit_predict_smoke.py" \
  --prepared "${PREPARED}/prepared" \
  --component-path "${COMPONENT}" \
  --helper-path "${ROOT}/scripts/multivi_fit_predict_smoke.py" \
  --output "${STAGE}/predictions" \
  --seed "${SCBUTTERFLY_SEED}" --epochs 30 --fold "${SCBUTTERFLY_FOLD}" \
  > "${STAGE}/validation/fit_predict.stdout.json" \
  2> "${STAGE}/validation/fit_predict.stderr.txt"
PYTHONPATH="${SITE_PACKAGES}" python -m pip freeze --all \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/scbutterfly_fit_predict_smoke.py" \
  "${ROOT}/scripts/multivi_fit_predict_smoke.py" \
  "${ROOT}/tests/unit/test_scbutterfly_fit_predict_smoke.py" \
  "${ROOT}/slurm/scbutterfly_fit_predict_fold.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${STAGE}" "${SCBUTTERFLY_FOLD}" "${SCBUTTERFLY_SEED}" <<'PY'
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
bundle_path = root / f"predictions/fold_{fold}/prediction_bundle.json"
bundle = json.loads(bundle_path.read_text())
if (
    receipt["status"] != "pass"
    or receipt["model_id"] != "scbutterfly_b_safe"
    or receipt["exact_distribution"] != "scButterfly 0.0.9"
    or receipt["outer_unit"] != "donor"
    or receipt["folds_requested"] != [fold]
    or set(receipt["folds"]) != {str(fold)}
    or receipt["query_atac_state"] != "structurally_missing"
    or receipt["held_atac_input_exposed"]
    or receipt["high_level_combined_preprocessing_executed"]
    or receipt["inverse_tfidf_executed"]
    or receipt["test_outcomes_read"]
    or receipt["champion_claim_allowed"]
    or not record.get("strict_state_resume")
    or not record.get("repeated_query_prediction_bit_identical")
    or record.get("held_atac_input_exposed")
    or record.get("inverse_tfidf_executed")
    or bundle.get("metadata", {}).get("seed") != seed + fold
):
    raise SystemExit("scButterfly fold receipt differs")
PredictionBundle.load_json(bundle_path)
freeze_tree(root, {
    "artifact_class": "scbutterfly_b_safe_fit_predict_fold",
    "model_id": "scbutterfly_b_safe",
    "dataset_id": "gse296875",
    "fold": fold,
    "base_seed": seed,
    "fold_seed": seed + fold,
    "prepared_artifacts_sha256": "6c10206033545fe7ed1ecbd9f7d64e107e2b3742de06db55a838bbdb46a1cdaa",
    "runtime_artifacts_sha256": "b495d6a28bff02387d4c1c56993df29bbac1cdc28119e8ec9bc4cff621965014",
    "held_atac_input_exposed": False,
    "test_outcomes_read": False,
    "champion_claim_allowed": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scButterfly donor-held fold: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
