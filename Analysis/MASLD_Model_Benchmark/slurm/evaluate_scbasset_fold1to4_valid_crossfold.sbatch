#!/usr/bin/env bash
#SBATCH --job-name=model-io-eval-032
#SBATCH --partition=io
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-io-eval-032-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-io-eval-032-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE_LOCK=${ROOT}/executions/scbasset-crossfold-evaluator-source-21070472
SOURCE_LOCK_SHA=be05e2b17d22de121880ca91b3b9750b0800f049731c05d1e243d7ad0441e9ad
PREDICTION_LOCK=${ROOT}/executions/scbasset-fold1to4-prediction-lock-21070492
PREDICTION_LOCK_SHA=9686a37c7bbeeac898cc2c4e1162c7033fcf2e0d09088abc9bb6c6c97bef3f07
EVALUATOR=${ROOT}/scripts/evaluate_scbasset_valid_crossfold.py
EVALUATOR_SHA=eec2e1f4e3f8d72a915af3f24a2959297b598a59aa2fdd0212c0019639553168
BASE_EVALUATOR=${ROOT}/scripts/evaluate_scbasset_valid.py
BASE_EVALUATOR_SHA=61018950b7772b306935f38ef0fa074ef5552aafa8fa71ce03bfa4c162676b3c
TEST=${ROOT}/tests/unit/test_evaluate_scbasset_valid_crossfold.py
TEST_SHA=3efad63b5a941ba963ccf5464014375a7e296d0a0d8d1cda0fdc8ff48bd10099
BIGWIGS=${ROOT}/executions/gse296875-deduplicated-tn5-bigwigs-21063992
BIGWIGS_SHA=01aacc2bb59898d6f7c133081b00c535047e271a2fed5f15761c2631b7545f1f
SPLIT=${ROOT}/executions/sequence-split-contract-v2-21063827
SPLIT_SHA=00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f
BASE=${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230
IMAGE=${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif
CORE=${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507
TARGET=${ROOT}/executions/scbasset-fold1to4-valid-crossfold-evaluation-${SLURM_JOB_ID}
STAGE=${TARGET}.staging

preserve_failed_stage() {
  local failed=${TARGET}.failed
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed scBasset cross-fold evaluation: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite scBasset cross-fold evaluation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
printf '%s  %s\n' \
  "${SOURCE_LOCK_SHA}" "${SOURCE_LOCK}/ARTIFACTS.json" \
  "${PREDICTION_LOCK_SHA}" "${PREDICTION_LOCK}/ARTIFACTS.json" \
  "${EVALUATOR_SHA}" "${EVALUATOR}" \
  "${BASE_EVALUATOR_SHA}" "${BASE_EVALUATOR}" \
  "${TEST_SHA}" "${TEST}" \
  b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 "${BASE}/ARTIFACTS.json" \
  96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa "${IMAGE}" \
  61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 "${CORE}/ARTIFACTS.json" \
  | sha256sum --check --strict

# Full prediction-side verification completes before outcome authorities are
# resolved or verified. The prediction record itself contains no outcome paths.
PYTHONPATH=${ROOT}/src python - "${SOURCE_LOCK}" "${PREDICTION_LOCK}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

source_lock = Path(sys.argv[1]).resolve(strict=True)
prediction_lock = Path(sys.argv[2]).resolve(strict=True)
verify_frozen_tree(source_lock)
verify_frozen_tree(prediction_lock)
payload = json.loads((prediction_lock / "prediction_lock.json").read_text(encoding="utf-8"))
if (
    payload["status"] != "pass"
    or payload["folds"] != [1, 2, 3, 4]
    or payload["outcomes_read"]
    or payload["outcome_paths_recorded"]
    or len(payload["records"]) != 4
):
    raise SystemExit("prediction lock contract differs before outcome access")
for record in payload["records"]:
    for key in ("input_root", "model_root", "prediction_root"):
        verify_frozen_tree(Path(record[key]))
PY

# Outcome authorities are included only after the prediction-side command
# above exits successfully.
printf '%s  %s\n' \
  "${BIGWIGS_SHA}" "${BIGWIGS}/ARTIFACTS.json" \
  "${SPLIT_SHA}" "${SPLIT}/ARTIFACTS.json" \
  | sha256sum --check --strict
PYTHONPATH=${ROOT}/src python - "${BIGWIGS}" "${SPLIT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
PY

module purge
module load singularity/4.2.2
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 \
  --env MKL_NUM_THREADS=1 \
  "${IMAGE}" python -m unittest -v \
  tests.unit.test_evaluate_scbasset_valid_crossfold \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 \
  --env MKL_NUM_THREADS=1 \
  "${IMAGE}" python "${EVALUATOR}" \
    --prediction-lock "${PREDICTION_LOCK}" \
    --prediction-lock-artifacts-sha256 "${PREDICTION_LOCK_SHA}" \
    --donor-bigwigs "${BIGWIGS}" \
    --donor-bigwigs-artifacts-sha256 "${BIGWIGS_SHA}" \
    --split-root "${SPLIT}" \
    --split-artifacts-sha256 "${SPLIT_SHA}" \
    --output "${STAGE}/evaluation" \
    --workers 8 \
    > "${STAGE}/validation/evaluation.stdout.json"
sha256sum "${EVALUATOR}" "${BASE_EVALUATOR}" "${TEST}" "$0" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}/src python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "evaluation/evaluation.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass"
    or receipt["evaluation_role"] != "valid"
    or receipt["outer_folds"] != [1, 2, 3, 4]
    or receipt["missing_outer_fold"] != 0
    or receipt["n_donors"] != 30
    or receipt["n_lineages"] != 5
    or receipt["n_regions_per_fold"] != 16000
    or not receipt["prediction_lock_verified_before_outcomes"]
    or receipt["held_valid_atac_exposed_to_model"]
    or receipt["test_regions_read"]
    or receipt["test_atac_read"]
    or receipt["test_predictions_read"]
    or receipt["champion_claim_allowed"]
):
    raise SystemExit("scBasset cross-fold evaluation contract differs")
print(freeze_tree(root, {
    "artifact_class": "scbasset_crossfold_valid_development_evaluation",
    "model_id": "scbasset",
    "dataset_id": "gse296875",
    "seed": 20260824,
    "outer_folds": [1, 2, 3, 4],
    "missing_outer_fold": 0,
    "n_donors": 30,
    "source_lock_artifacts_sha256": "be05e2b17d22de121880ca91b3b9750b0800f049731c05d1e243d7ad0441e9ad",
    "prediction_lock_artifacts_sha256": receipt["prediction_lock_artifacts_sha256"],
    "held_valid_outcomes_read": True,
    "test_outcomes_read": False,
    "champion_claim_allowed": False,
    "status": "pass",
}))
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scBasset cross-fold validation evaluation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
