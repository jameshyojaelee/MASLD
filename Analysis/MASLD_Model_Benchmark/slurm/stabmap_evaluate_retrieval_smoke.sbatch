#!/usr/bin/env bash
#SBATCH --job-name=masld-stabmap-eval
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
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

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
FOLD0="${ROOT}/executions/stabmap-retrieval-fold0-21066517"
FOLD1="${ROOT}/executions/stabmap-retrieval-fold1-21066516"
FOLD2="${ROOT}/executions/stabmap-retrieval-fold2-21066515"
FOLD3="${ROOT}/executions/stabmap-retrieval-fold3-21066518"
FOLD4="${ROOT}/executions/stabmap-retrieval-fold4-21066490"
BASELINE="${ROOT}/executions/scglue-retrieval-fit-predict-21065540"
BASELINE_ARTIFACTS_SHA256=6a37d1c78d1bf6d798e68d8cb6d4e683befbedc4d1a21eea6ec08fb9c82656ad
PREPARED="${ROOT}/executions/scglue-retrieval-prepared-21065511"
PREPARED_ARTIFACTS_SHA256=c79930c961cbbd44d9d198579f8b9443b6758e28cc2f4458be48d884d62a3354
RUNTIME="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
RUNTIME_ARTIFACTS_SHA256=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
SITE_PACKAGES="${RUNTIME}/site-packages"
TARGET="${ROOT}/executions/stabmap-retrieval-evaluation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed StabMap retrieval evaluation: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite StabMap retrieval evaluation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for record in \
  "e9dd7f9cf4ec533d62f5b742a4f10241506216b9828f8f9241a93f1ea85476bf ${FOLD0}/ARTIFACTS.json" \
  "0219c04403b2b87ad27202eb04dcf84385e346ee6c63c10144d78c6c2d45de82 ${FOLD1}/ARTIFACTS.json" \
  "2767de0cb14867017750dee4904dd74142ee75ae5d80126eaa5450120fcfc11c ${FOLD2}/ARTIFACTS.json" \
  "96a82dfda4f805f5594e1e16ae3e965dd381b1670bc9c1c363b5764d6fdb8d5c ${FOLD3}/ARTIFACTS.json" \
  "1a9f775d9e2fc51e55e3b5864b56b7b0ee7c252675582e0796105b54e969b182 ${FOLD4}/ARTIFACTS.json" \
  "${BASELINE_ARTIFACTS_SHA256} ${BASELINE}/ARTIFACTS.json" \
  "${PREPARED_ARTIFACTS_SHA256} ${PREPARED}/ARTIFACTS.json" \
  "${RUNTIME_ARTIFACTS_SHA256} ${RUNTIME}/ARTIFACTS.json"; do
  printf '%s\n' "${record}" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${FOLD0}" "${FOLD1}" "${FOLD2}" "${FOLD3}" "${FOLD4}" \
  "${BASELINE}" "${PREPARED}" "${RUNTIME}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for path in sys.argv[1:]:
    verify_frozen_tree(Path(path))
PY
module purge
module load python/3.12.3-GCCcore-13.3.0
PYTHONPATH="${SITE_PACKAGES}:${ROOT}" python -m unittest -v \
  tests.unit.test_evaluate_scglue_retrieval_smoke \
  tests.unit.test_evaluate_stabmap_retrieval_smoke \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${SITE_PACKAGES}:${ROOT}/scripts" python \
  "${ROOT}/scripts/evaluate_stabmap_retrieval_smoke.py" \
  --fold-root "0=${FOLD0}" --fold-root "1=${FOLD1}" \
  --fold-root "2=${FOLD2}" --fold-root "3=${FOLD3}" \
  --fold-root "4=${FOLD4}" \
  --baseline-predictions "${BASELINE}/predictions" \
  --evaluator-only "${PREPARED}/prepared/evaluator_only" \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/validation/evaluation.stdout.json"
sha256sum \
  "${ROOT}/scripts/evaluate_stabmap_retrieval_smoke.py" \
  "${ROOT}/scripts/evaluate_scglue_retrieval_smoke.py" \
  "${ROOT}/tests/unit/test_evaluate_stabmap_retrieval_smoke.py" \
  "${ROOT}/tests/unit/test_evaluate_scglue_retrieval_smoke.py" \
  "${ROOT}/slurm/stabmap_evaluate_retrieval_smoke.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "evaluation/evaluation.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["model_roster"] != ["stabmap", "linear_cca"]
    or receipt["task_id"] != "same_nucleus_cross_modal_retrieval"
    or receipt["outer_unit"] != "donor"
    or receipt["n_donors"] != 39
    or receipt["hidden_pair_map_available_to_model"]
    or receipt["cells_used_as_biological_replicates"]
    or receipt["test_outcomes_read"]
    or receipt["rna_to_atac_prediction_claim"]
    or receipt["sealed_rna_conditioned_atac_eligible"]
    or receipt["champion_claim_allowed"]
):
    raise SystemExit("StabMap retrieval evaluation contract differs")
freeze_tree(root, {
    "artifact_class": "stabmap_retrieval_development_evaluation",
    "model_id": "stabmap",
    "dataset_id": "gse296875",
    "fold_artifact_sha256": {
        "0": "e9dd7f9cf4ec533d62f5b742a4f10241506216b9828f8f9241a93f1ea85476bf",
        "1": "0219c04403b2b87ad27202eb04dcf84385e346ee6c63c10144d78c6c2d45de82",
        "2": "2767de0cb14867017750dee4904dd74142ee75ae5d80126eaa5450120fcfc11c",
        "3": "96a82dfda4f805f5594e1e16ae3e965dd381b1670bc9c1c363b5764d6fdb8d5c",
        "4": "1a9f775d9e2fc51e55e3b5864b56b7b0ee7c252675582e0796105b54e969b182",
    },
    "baseline_artifacts_sha256": "6a37d1c78d1bf6d798e68d8cb6d4e683befbedc4d1a21eea6ec08fb9c82656ad",
    "hidden_pair_map_available_to_model": False,
    "cells_used_as_biological_replicates": False,
    "test_outcomes_read": False,
    "sealed_rna_conditioned_atac_eligible": False,
    "champion_claim_allowed": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'StabMap retrieval development evaluation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
