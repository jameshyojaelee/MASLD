#!/usr/bin/env bash
#SBATCH --job-name=model-data-053
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
RAW=${ROOT}/executions/geneformer-v2-316m-frozen-screen-extraction-21082667
RAW_SHA256=06cf29415a5e6696e235c2f9a5ffc4edd8127043327f6c36e237e8d930390411
SPLIT=${ROOT}/executions/model-data-052-21077595/split
SPLIT_SHA256=10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680
CHECKPOINT_REGISTRY=${ROOT}/config/artifacts/models/geneformer/checkpoints.json
CHECKPOINT_REGISTRY_SHA256=11f311a5078fc4c8f9a98e55fea6b4ccf278b2b01503da8086ac3f90cdc73455
EXPOSURE_AUDIT=${ROOT}/config/artifacts/models/geneformer/exposure_audit.json
EXPOSURE_AUDIT_SHA256=6aef01d60a68689af920a5c571b0979bdd7364d39d35e59cc1a2ba3e242032d3
TARGET=${ROOT}/executions/model-data-053-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite common embedding preparation artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${RAW_SHA256}" "${RAW}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SPLIT_SHA256}" "${SPLIT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${CHECKPOINT_REGISTRY_SHA256}" "${CHECKPOINT_REGISTRY}" \
  | sha256sum --check --strict
printf '%s  %s\n' "${EXPOSURE_AUDIT_SHA256}" "${EXPOSURE_AUDIT}" \
  | sha256sum --check --strict

cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
"${PYTHON}" -m unittest -v \
  tests.unit.test_prepare_common_cell_embeddings_50000 \
  tests.unit.test_audit_evaluate_common_cell_heads_study_50000 \
  tests.unit.test_fit_predict_common_cell_heads_study_50000 \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${PYTHON}" scripts/prepare_common_cell_embeddings_50000.py \
  --raw-root "${RAW}" \
  --raw-embeddings-relative extraction/embeddings.npy \
  --raw-row-order-relative extraction/embedding_row_order.txt \
  --raw-receipt-relative extraction/extraction_receipt.json \
  --split-root "${SPLIT}" \
  --checkpoint-registry "${CHECKPOINT_REGISTRY}" \
  --exposure-audit "${EXPOSURE_AUDIT}" \
  --output "${STAGE}/embeddings" \
  --expected-raw-artifacts-sha256 "${RAW_SHA256}" \
  --expected-split-artifacts-sha256 "${SPLIT_SHA256}" \
  --expected-checkpoint-registry-sha256 "${CHECKPOINT_REGISTRY_SHA256}" \
  --expected-exposure-audit-sha256 "${EXPOSURE_AUDIT_SHA256}" \
  --expected-model-id geneformer_v2_316m \
  --expected-checkpoint-artifact-path Geneformer-V2-316M/model.safetensors \
  --expected-checkpoint-sha256 965ceccea81953d362081ef3843560a0e4fef88d396c28017881f1e94b1246f3 \
  --expected-exposure-status target_label_unexposed \
  --expected-sealed-champion-eligible true \
  --expected-rows 50000 \
  --expected-width 1152 \
  --expected-donors 102 \
  --expected-studies 7 \
  > "${STAGE}/validation/preparation_stdout.json" \
  2> "${STAGE}/validation/preparation_stderr.txt"
"${PYTHON}" - "${STAGE}/embeddings" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree
from scripts.fit_predict_common_cell_heads_study_50000 import derive_embedding_identity

root = Path(sys.argv[1]).resolve(strict=True)
manifest = verify_frozen_tree(root)
identity, path = derive_embedding_identity(
    root, Path("embeddings/common_embeddings.npz")
)
receipt = json.loads((root / "embeddings/receipt.json").read_text())
if (
    manifest["metadata"].get("artifact_class")
    != "cell_foundation_common_embeddings"
    or identity
    != {
        "model_id": "geneformer_v2_316m",
        "representation_id": "common",
        "checkpoint_sha256": "965ceccea81953d362081ef3843560a0e4fef88d396c28017881f1e94b1246f3",
        "exposure_status": "target_label_unexposed",
        "embedding_width": 1152,
        "embedding_relative": "embeddings/common_embeddings.npz",
        "sealed_champion_eligible": True,
    }
    or path.stat().st_size < 230_000_000
    or receipt.get("rows") != 50_000
    or receipt.get("source_label_table_read") is not False
    or receipt.get("evaluation_label_columns_read") != []
    or receipt.get("downstream_head_fit") is not False
    or receipt.get("histology_read") is not False
    or receipt.get("sealed_outcomes_read") is not False
):
    raise SystemExit("Geneformer common embedding contract differs")
PY
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/prepare_common_cell_embeddings_50000.py \
  scripts/audit_aggregate_common_cell_heads_study_50000.py \
  scripts/score_common_cell_heads_study_50000.py \
  scripts/fit_predict_common_cell_heads_study_50000.py \
  tests/unit/test_prepare_common_cell_embeddings_50000.py \
  tests/unit/test_audit_evaluate_common_cell_heads_study_50000.py \
  tests/unit/test_fit_predict_common_cell_heads_study_50000.py \
  slurm/prepare_geneformer_v2_316m_common_50000_cpu.sbatch \
  slurm/run_common_cell_head_bundle_50000.sbatch \
  slurm/audit_aggregate_common_cell_heads_study_50000.sbatch \
  slurm/score_common_cell_heads_study_50000.sbatch \
  "${RAW}/ARTIFACTS.json" \
  "${SPLIT}/ARTIFACTS.json" \
  "${CHECKPOINT_REGISTRY}" \
  "${EXPOSURE_AUDIT}" \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
manifest = verify_frozen_tree(root / "embeddings")
freeze_tree(
    root,
    {
        "artifact_class": "geneformer_v2_316m_common_embedding_preparation_campaign",
        "model_id": "geneformer_v2_316m",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "embedding_artifacts_sha256": __import__("hashlib").sha256(
            (root / "embeddings/ARTIFACTS.json").read_bytes()
        ).hexdigest(),
        "evaluation_labels_read": False,
        "downstream_head_fit": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/embeddings/ARTIFACTS.json"
