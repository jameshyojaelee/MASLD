#!/usr/bin/env bash
#SBATCH --job-name=masld-midas-aggregate
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

: "${MIDAS_SEED:?MIDAS_SEED is required}"
: "${MIDAS_FOLD0:?MIDAS_FOLD0 is required}"
: "${MIDAS_FOLD1:?MIDAS_FOLD1 is required}"
: "${MIDAS_FOLD2:?MIDAS_FOLD2 is required}"
: "${MIDAS_FOLD3:?MIDAS_FOLD3 is required}"
: "${MIDAS_FOLD4:?MIDAS_FOLD4 is required}"
: "${MIDAS_FOLD0_SHA:?MIDAS_FOLD0_SHA is required}"
: "${MIDAS_FOLD1_SHA:?MIDAS_FOLD1_SHA is required}"
: "${MIDAS_FOLD2_SHA:?MIDAS_FOLD2_SHA is required}"
: "${MIDAS_FOLD3_SHA:?MIDAS_FOLD3_SHA is required}"
: "${MIDAS_FOLD4_SHA:?MIDAS_FOLD4_SHA is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/midas-inductive-seed${MIDAS_SEED}-aggregated-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed MIDAS aggregation: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite MIDAS aggregation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for record in \
  "${MIDAS_FOLD0_SHA} ${MIDAS_FOLD0}/ARTIFACTS.json" \
  "${MIDAS_FOLD1_SHA} ${MIDAS_FOLD1}/ARTIFACTS.json" \
  "${MIDAS_FOLD2_SHA} ${MIDAS_FOLD2}/ARTIFACTS.json" \
  "${MIDAS_FOLD3_SHA} ${MIDAS_FOLD3}/ARTIFACTS.json" \
  "${MIDAS_FOLD4_SHA} ${MIDAS_FOLD4}/ARTIFACTS.json"; do
  printf '%s\n' "${record}" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${MIDAS_FOLD0}" "${MIDAS_FOLD1}" "${MIDAS_FOLD2}" \
  "${MIDAS_FOLD3}" "${MIDAS_FOLD4}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for path in sys.argv[1:]:
    verify_frozen_tree(Path(path))
PY
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.unit.test_aggregate_midas_fold_predictions \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}/scripts" python \
  "${ROOT}/scripts/aggregate_midas_fold_predictions.py" \
  --fold-root "0=${MIDAS_FOLD0}" --fold-root "1=${MIDAS_FOLD1}" \
  --fold-root "2=${MIDAS_FOLD2}" --fold-root "3=${MIDAS_FOLD3}" \
  --fold-root "4=${MIDAS_FOLD4}" \
  --expected-seed "${MIDAS_SEED}" --output "${STAGE}/predictions" \
  > "${STAGE}/validation/aggregation.stdout.json" \
  2> "${STAGE}/validation/aggregation.stderr.txt"
sha256sum \
  "${ROOT}/scripts/aggregate_midas_fold_predictions.py" \
  "${ROOT}/tests/unit/test_aggregate_midas_fold_predictions.py" \
  "${ROOT}/slurm/midas_aggregate_folds.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - \
  "${STAGE}" "${MIDAS_SEED}" \
  "${MIDAS_FOLD0_SHA}" "${MIDAS_FOLD1_SHA}" "${MIDAS_FOLD2_SHA}" \
  "${MIDAS_FOLD3_SHA}" "${MIDAS_FOLD4_SHA}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench import PredictionBundle
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
seed = int(sys.argv[2])
fold_sha = {str(fold): value for fold, value in enumerate(sys.argv[3:8])}
receipt = json.loads((root / "predictions/receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["model_id"] != "midas_inductive"
    or receipt["fold_count"] != 5
    or set(receipt["folds"]) != set(map(str, range(5)))
    or receipt["base_seed"] != seed
    or not receipt["prediction_frozen_before_outcomes"]
    or receipt["query_rna_used_for_training"]
    or receipt["held_atac_input_exposed"]
    or receipt["outcomes_read"]
    or receipt["test_outcomes_read"]
    or receipt["sealed_inference_eligible"]
    or receipt["champion_claim_allowed"]
):
    raise SystemExit("MIDAS aggregation receipt differs")
for fold in range(5):
    PredictionBundle.load_json(root / f"predictions/fold_{fold}/prediction_bundle.json")
freeze_tree(root, {
    "artifact_class": "midas_inductive_aggregated_predictions",
    "model_id": "midas_inductive",
    "dataset_id": "gse296875",
    "seed": seed,
    "fold_artifacts_sha256": fold_sha,
    "prediction_frozen_before_outcomes": True,
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
printf 'MIDAS aggregated predictions: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
