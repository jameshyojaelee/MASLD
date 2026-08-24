#!/usr/bin/env bash
#SBATCH --job-name=masld-seurat-wnn-fit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
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
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PREPARED="${ROOT}/executions/seurat-wnn-bridge-prepared-21065772"
PREPARED_ARTIFACTS_SHA256=099bea7c9a3728da1eb0eaf2938f1a8f26fa9de6f6c087532aaadc4df9ee3500
RUNTIME="${ROOT}/executions/environments/seurat-5.5.1-portable-21064881"
RUNTIME_ARTIFACTS_SHA256=67f5689231eac5f66cb82b830c32217943071b59b74180fdc1031978b3478eff
TARGET="${ROOT}/executions/seurat-wnn-bridge-fit-predict-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed Seurat WNN bridge fit: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite Seurat WNN bridge fit: %s\n' "${TARGET}" >&2
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
module purge
module load R/4.5.2-mba
R_LIBS="${RUNTIME}/library" Rscript --vanilla \
  "${ROOT}/scripts/seurat_wnn_bridge_fit_predict_smoke.R" --self-test \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
R_LIBS="${RUNTIME}/library" Rscript --vanilla \
  "${ROOT}/scripts/seurat_wnn_bridge_fit_predict_smoke.R" \
  "${PREPARED}/prepared" "${STAGE}/predictions" \
  > "${STAGE}/validation/fit_predict.stdout.json" \
  2> "${STAGE}/validation/fit_predict.stderr.txt"
R_LIBS="${RUNTIME}/library" Rscript --vanilla - <<'RS' \
  > "${STAGE}/validation/sessionInfo.txt"
sessionInfo()
RS
sha256sum \
  "${ROOT}/scripts/seurat_wnn_bridge_fit_predict_smoke.R" \
  "${ROOT}/slurm/seurat_wnn_bridge_fit_predict_smoke.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "predictions/receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["model_id"] != "seurat_wnn_bridge"
    or receipt["exact_package"] != "Seurat 5.5.1"
    or receipt["pairing_topology"] != "same_nucleus"
    or receipt["outer_unit"] != "donor"
    or set(receipt["folds"]) != set(map(str, range(5)))
    or not receipt["training_pairing_available"]
    or not receipt["query_modality_ids_disjoint"]
    or receipt["hidden_pair_map_read"]
    or receipt["retrieval_metrics_calculated"]
    or receipt["cells_used_as_biological_replicates"]
    or receipt["outcomes_read"]
    or receipt["test_outcomes_read"]
    or receipt["rna_to_atac_prediction_executed"]
    or receipt["sealed_rna_conditioned_atac_eligible"]
    or receipt["champion_claim_allowed"]
):
    raise SystemExit("Seurat WNN bridge fit contract differs")
for fold, value in receipt["folds"].items():
    if (
        not value["strict_state_resume"]
        or not value["repeated_query_mapping_bit_identical"]
        or value["hidden_pair_map_read"]
        or value["retrieval_metrics_calculated"]
    ):
        raise SystemExit(f"Seurat WNN bridge fold contract differs: {fold}")
freeze_tree(root, {
    "artifact_class": "seurat_wnn_bridge_fit_predict",
    "model_id": "seurat_wnn_bridge",
    "exact_package": "Seurat 5.5.1",
    "dataset_id": "gse296875",
    "prepared_artifacts_sha256": "099bea7c9a3728da1eb0eaf2938f1a8f26fa9de6f6c087532aaadc4df9ee3500",
    "runtime_artifacts_sha256": "67f5689231eac5f66cb82b830c32217943071b59b74180fdc1031978b3478eff",
    "hidden_pair_map_read": False,
    "retrieval_metrics_calculated": False,
    "cells_used_as_biological_replicates": False,
    "outcomes_read": False,
    "test_outcomes_read": False,
    "rna_to_atac_prediction_executed": False,
    "sealed_rna_conditioned_atac_eligible": False,
    "champion_claim_allowed": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Seurat WNN bridge fit/predict: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
