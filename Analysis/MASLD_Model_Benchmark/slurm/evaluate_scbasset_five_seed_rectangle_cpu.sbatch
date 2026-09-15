#!/usr/bin/env bash
#SBATCH --job-name=model-eval
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-eval-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-eval-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

: "${COMPLETION_ROOT:?completion root must be exported by the exact-coverage gate}"
: "${COMPLETION_SHA:?completion ARTIFACTS hash must be exported by the exact-coverage gate}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONTRACT=${ROOT}/config/evaluation/scbasset_five_seed_rectangular_evaluator_20260825.json
BASE_EVALUATOR=${ROOT}/scripts/evaluate_scbasset_valid.py
EVALUATOR=${ROOT}/scripts/evaluate_scbasset_five_seed_rectangle.py
ACCEPTANCE=${ROOT}/executions/sequence-regulatory-completion-gate-tests-21100975
ENVIRONMENT=${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230
IMAGE=${ENVIRONMENT}/tensorflow-tensorflow-2.8.2-gpu.sif
CORE=${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507
TARGET=${ROOT}/executions/scbasset-five-seed-evaluation-${SLURM_JOB_ID}
STAGE=${TARGET}.staging

preserve_failed_stage() {
  local failed=${TARGET}.failed
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite scBasset evaluation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

cd "${ROOT}"
printf '%s  %s\n' \
  69a7f5a504476482b2e32a9486a9f37b43b4e3cb3cf297422026303be03b39da "${CONTRACT}" \
  61018950b7772b306935f38ef0fa074ef5552aafa8fa71ce03bfa4c162676b3c "${BASE_EVALUATOR}" \
  1c314b6aee1ec0d8043fd8a1221c8ee69b5781bb240acbf78d4b8a30feb353e6 "${EVALUATOR}" \
  7c230d85f1caafff06a2fc1be7bf79a5e5dc3857b2f93d1833f17893e8f37027 tests/unit/test_evaluate_scbasset_five_seed_rectangle.py \
  4632e8d6debea40b832365c31fffd301bcd9896dd765f2722f3b2332ff1e4d29 "${ACCEPTANCE}/ARTIFACTS.json" \
  b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 "${ENVIRONMENT}/ARTIFACTS.json" \
  96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa "${IMAGE}" \
  61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 "${CORE}/ARTIFACTS.json" \
  "${COMPLETION_SHA}" "${COMPLETION_ROOT}/ARTIFACTS.json" \
  | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}/src:${ROOT} python - "${ACCEPTANCE}" "${COMPLETION_ROOT}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import verify_frozen_tree

acceptance, completion = map(lambda value: Path(value).resolve(strict=True), sys.argv[1:])
verify_frozen_tree(acceptance)
verify_frozen_tree(completion)
receipt = json.loads(
    (completion / "completion/completion_lock.json").read_text(encoding="utf-8")
)
if (
    receipt["status"] != "ready_both_rectangles_locked"
    or not receipt["evaluation_launch_authorized"]
    or receipt["partial_ranking_authorized"]
    or receipt["prediction_values_read"]
    or receipt["observed_outcomes_read"]
):
    raise SystemExit("exact-coverage completion lock differs")
PY

module purge
module load singularity/4.2.2
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${ROOT}/src:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 \
  --env MKL_NUM_THREADS=1 \
  "${IMAGE}" python "${EVALUATOR}" \
    --root "${ROOT}" \
    --evaluator-contract "${CONTRACT}" \
    --completion-root "${COMPLETION_ROOT}" \
    --completion-artifacts-sha256 "${COMPLETION_SHA}" \
    --output "${STAGE}/evaluation" \
    --workers 8 \
    --bootstrap-resamples 10000 \
    --bootstrap-seed 17 \
    > "${STAGE}/validation/evaluation.stdout.json" \
    2> "${STAGE}/validation/evaluation.stderr.txt"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${ROOT}/src:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${CONTRACT}" "${BASE_EVALUATOR}" "${EVALUATOR}" \
  tests/unit/test_evaluate_scbasset_five_seed_rectangle.py \
  "$0" > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}/src:${ROOT} python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
result = json.loads((root / "evaluation/evaluation.json").read_text(encoding="utf-8"))
if (
    result["status"] != "pass"
    or result["fits"] != 25
    or result["prediction_views"] != 25
    or result["n_donors"] != 39
    or result["partial_rectangle_used"]
    or result["cross_task_family_ranking_performed"]
    or result["held_valid_atac_exposed_to_model"]
    or result["test_outcomes_read"]
    or result["sealed_data_read"]
    or result["rna_conditioned_claim_allowed"]
    or result["champion_claim_allowed"]
    or result["universal_claim_allowed"]
):
    raise SystemExit("scBasset evaluation terminal contract differs")
print(freeze_tree(root, {
    "artifact_class": "scbasset_five_seed_development_evaluation",
    "dataset_id": "gse296875",
    "status": "pass",
    "fits": 25,
    "prediction_views": 25,
    "n_donors": 39,
    "completion_artifacts_sha256": result["completion_artifacts_sha256"],
    "partial_rectangle_used": False,
    "cross_task_family_ranking_performed": False,
    "rna_conditioned_claim_allowed": False,
    "test_outcomes_read": False,
    "sealed_data_read": False,
    "champion_claim_allowed": False,
    "universal_claim_allowed": False,
}))
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scBasset five-seed evaluation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
