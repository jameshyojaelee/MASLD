#!/usr/bin/env bash
#SBATCH --job-name=masld-stabmap-fold
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export R_ENVIRON_USER=/dev/null
export R_PROFILE_USER=/dev/null
export R_LIBS_USER=/nonexistent/masld_bench_r_user_library
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

if [[ ! "${FOLD:-}" =~ ^[0-4]$ ]]; then
  printf 'FOLD must be one of 0,1,2,3,4\n' >&2
  exit 2
fi
if [[ ! "${STABMAP_SEED_BASE:-4817}" =~ ^[0-9]+$ ]]; then
  printf 'STABMAP_SEED_BASE must be a nonnegative integer\n' >&2
  exit 2
fi
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PREPARED="${ROOT}/executions/seurat-wnn-bridge-prepared-21065772"
PREPARED_ARTIFACTS_SHA256=099bea7c9a3728da1eb0eaf2938f1a8f26fa9de6f6c087532aaadc4df9ee3500
RUNTIME="${ROOT}/executions/environments/stabmap-v1.6.0-sparse-portable-21066444"
RUNTIME_ARTIFACTS_SHA256=ec6963d24e95378073e5e27c2305ce854ed11fc2faeba8efc8638bc8dca71cf5
TARGET="${ROOT}/executions/stabmap-retrieval-seed${STABMAP_SEED_BASE}-fold${FOLD}-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed StabMap fold: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite StabMap fold: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for record in \
  "${PREPARED_ARTIFACTS_SHA256} ${PREPARED}/ARTIFACTS.json" \
  "${RUNTIME_ARTIFACTS_SHA256} ${RUNTIME}/ARTIFACTS.json"; do
  printf '%s\n' "${record}" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${PREPARED}" "${RUNTIME}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for path in sys.argv[1:]:
    verify_frozen_tree(Path(path))
PY
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.unit.test_stabmap_fit_predict_retrieval_fold \
  > "${STAGE}/validation/unit.stdout.txt" \
  2> "${STAGE}/validation/unit.stderr.txt"
module purge
module load R/4.5.2-mba
R_LIBS="${RUNTIME}/library" Rscript --vanilla \
  "${ROOT}/scripts/stabmap_fit_predict_retrieval_fold.R" --self-test \
  > "${STAGE}/validation/self_test.stdout.txt" \
  2> "${STAGE}/validation/self_test.stderr.txt"
R_LIBS="${RUNTIME}/library" Rscript --vanilla \
  "${ROOT}/scripts/stabmap_fit_predict_retrieval_fold.R" \
  "${PREPARED}/prepared/fold_${FOLD}" "${STAGE}/predictions/fold_${FOLD}/stabmap" \
  "${FOLD}" "${STABMAP_SEED_BASE}" \
  > "${STAGE}/validation/fit_predict.stdout.tsv" \
  2> "${STAGE}/validation/fit_predict.stderr.txt"
R_LIBS="${RUNTIME}/library" Rscript --vanilla - <<'RS' \
  > "${STAGE}/validation/sessionInfo.txt"
sessionInfo()
RS
sha256sum \
  "${ROOT}/scripts/stabmap_fit_predict_retrieval_fold.R" \
  "${ROOT}/tests/unit/test_stabmap_fit_predict_retrieval_fold.py" \
  "${ROOT}/slurm/stabmap_fit_predict_retrieval_fold.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${FOLD}" "${STABMAP_SEED_BASE}" <<'PY'
from pathlib import Path
import csv
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
fold = int(sys.argv[2])
base_seed = int(sys.argv[3])
receipt_path = root / f"predictions/fold_{fold}/stabmap/receipt.tsv"
with receipt_path.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
receipt = {row["key"]: row["value"] for row in rows}
expected = {
    "status": "pass",
    "model_id": "stabmap",
    "exact_package": "StabMap 1.6.0",
    "dataset_id": "gse296875",
    "fold": str(fold),
    "seed": str(base_seed + fold * 101),
    "embedding_dimensions": "16",
    "training_state_resume": "true",
    "query_order_invariant": "true",
    "hidden_pair_map_read": "false",
    "retrieval_metrics_calculated": "false",
    "outcomes_read": "false",
    "test_outcomes_read": "false",
}
if any(receipt.get(key) != value for key, value in expected.items()):
    raise SystemExit("StabMap fold receipt differs")
for modality in ("rna", "atac"):
    path = root / f"predictions/fold_{fold}/stabmap/query_{modality}_embeddings.tsv"
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        rows = list(reader)
    if len(header) != 17 or header[0] != "query_id" or len(rows) < 100:
        raise SystemExit(f"StabMap {modality} embedding axes differ")
freeze_tree(root, {
    "artifact_class": "stabmap_retrieval_fit_predict_fold",
    "model_id": "stabmap",
    "exact_package": "StabMap 1.6.0",
    "dataset_id": "gse296875",
    "fold": fold,
    "base_seed": base_seed,
    "fold_seed": base_seed + fold * 101,
    "prepared_artifacts_sha256": "099bea7c9a3728da1eb0eaf2938f1a8f26fa9de6f6c087532aaadc4df9ee3500",
    "runtime_artifacts_sha256": "ec6963d24e95378073e5e27c2305ce854ed11fc2faeba8efc8638bc8dca71cf5",
    "pairing_topology": "same_nucleus",
    "outer_unit": "donor",
    "terminal_queries_identity_blinded": True,
    "hidden_pair_map_read": False,
    "retrieval_metrics_calculated": False,
    "cells_used_as_biological_replicates": False,
    "outcomes_read": False,
    "test_outcomes_read": False,
    "sealed_rna_conditioned_atac_eligible": False,
    "champion_claim_allowed": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'StabMap fold %s: %s\n' "${FOLD}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
