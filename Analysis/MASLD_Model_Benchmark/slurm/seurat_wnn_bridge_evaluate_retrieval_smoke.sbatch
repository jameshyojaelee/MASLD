#!/usr/bin/env bash
#SBATCH --job-name=masld-seurat-wnn-eval
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
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
WNN="${ROOT}/executions/seurat-wnn-bridge-fit-predict-21065917"
WNN_ARTIFACTS_SHA256=f8fbfe21e08100ac86d70bb45566057b8325fe8a223df5c9f4ef41b6fbe06af3
BASELINE="${ROOT}/executions/scglue-retrieval-fit-predict-21065540"
BASELINE_ARTIFACTS_SHA256=6a37d1c78d1bf6d798e68d8cb6d4e683befbedc4d1a21eea6ec08fb9c82656ad
PREPARED="${ROOT}/executions/scglue-retrieval-prepared-21065511"
PREPARED_ARTIFACTS_SHA256=c79930c961cbbd44d9d198579f8b9443b6758e28cc2f4458be48d884d62a3354
RUNTIME="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
RUNTIME_ARTIFACTS_SHA256=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
SITE_PACKAGES="${RUNTIME}/site-packages"
TARGET="${ROOT}/executions/seurat-wnn-bridge-evaluation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed Seurat WNN bridge evaluation: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite Seurat WNN bridge evaluation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for record in \
  "${WNN_ARTIFACTS_SHA256} ${WNN}/ARTIFACTS.json" \
  "${BASELINE_ARTIFACTS_SHA256} ${BASELINE}/ARTIFACTS.json" \
  "${PREPARED_ARTIFACTS_SHA256} ${PREPARED}/ARTIFACTS.json" \
  "${RUNTIME_ARTIFACTS_SHA256} ${RUNTIME}/ARTIFACTS.json"; do
  printf '%s\n' "${record}" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${WNN}" "${BASELINE}" "${PREPARED}" "${RUNTIME}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for path in sys.argv[1:]:
    verify_frozen_tree(Path(path))
PY
module purge
module load python/3.12.3-GCCcore-13.3.0
PYTHONPATH="${SITE_PACKAGES}:${ROOT}/scripts:${ROOT}" python -m unittest -v \
  tests.unit.test_evaluate_seurat_wnn_bridge_retrieval_smoke \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${SITE_PACKAGES}:${ROOT}/scripts" python \
  "${ROOT}/scripts/evaluate_seurat_wnn_bridge_retrieval_smoke.py" \
  --wnn-predictions "${WNN}/predictions" \
  --baseline-predictions "${BASELINE}/predictions" \
  --evaluator-only "${PREPARED}/prepared/evaluator_only" \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/validation/evaluation.stdout.json"
sha256sum \
  "${ROOT}/scripts/evaluate_seurat_wnn_bridge_retrieval_smoke.py" \
  "${ROOT}/scripts/evaluate_scglue_retrieval_smoke.py" \
  "${ROOT}/tests/unit/test_evaluate_seurat_wnn_bridge_retrieval_smoke.py" \
  "${ROOT}/slurm/seurat_wnn_bridge_evaluate_retrieval_smoke.sbatch" \
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
    or receipt["task_id"] != "same_nucleus_cross_modal_retrieval"
    or receipt["model_roster"] != ["seurat_wnn_bridge", "linear_cca"]
    or receipt["outer_unit"] != "donor"
    or receipt["n_donors"] != 39
    or receipt["n_same_nucleus_pairs"] != 1000
    or receipt["hidden_pair_map_available_to_model"]
    or not receipt["query_modality_ids_disjoint"]
    or not receipt["query_atac_row_order_permuted"]
    or receipt["cells_used_as_biological_replicates"]
    or not receipt["donor_macro_aggregation"]
    or receipt["test_outcomes_read"]
    or receipt["rna_to_atac_prediction_claim"]
    or receipt["sealed_rna_conditioned_atac_eligible"]
    or receipt["champion_claim_allowed"]
):
    raise SystemExit("Seurat WNN bridge evaluation contract differs")
freeze_tree(root, {
    "artifact_class": "seurat_wnn_bridge_retrieval_development_evaluation",
    "model_roster": ["seurat_wnn_bridge", "linear_cca"],
    "dataset_id": "gse296875",
    "wnn_artifacts_sha256": "f8fbfe21e08100ac86d70bb45566057b8325fe8a223df5c9f4ef41b6fbe06af3",
    "baseline_artifacts_sha256": "6a37d1c78d1bf6d798e68d8cb6d4e683befbedc4d1a21eea6ec08fb9c82656ad",
    "prepared_artifacts_sha256": "c79930c961cbbd44d9d198579f8b9443b6758e28cc2f4458be48d884d62a3354",
    "hidden_pair_map_available_to_model": False,
    "cells_used_as_biological_replicates": False,
    "test_outcomes_read": False,
    "rna_to_atac_prediction_claim": False,
    "sealed_rna_conditioned_atac_eligible": False,
    "champion_claim_allowed": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Seurat WNN bridge retrieval evaluation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
