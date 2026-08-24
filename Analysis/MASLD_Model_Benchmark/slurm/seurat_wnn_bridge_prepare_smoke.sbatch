#!/usr/bin/env bash
#SBATCH --job-name=masld-seurat-wnn-prepare
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
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
SOURCE="${ROOT}/executions/scglue-retrieval-prepared-21065511"
SOURCE_ARTIFACTS_SHA256=c79930c961cbbd44d9d198579f8b9443b6758e28cc2f4458be48d884d62a3354
RUNTIME="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
RUNTIME_ARTIFACTS_SHA256=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
SITE_PACKAGES="${RUNTIME}/site-packages"
TARGET="${ROOT}/executions/seurat-wnn-bridge-prepared-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed Seurat WNN bridge preparation: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite Seurat WNN bridge preparation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
for record in \
  "${SOURCE_ARTIFACTS_SHA256} ${SOURCE}/ARTIFACTS.json" \
  "${RUNTIME_ARTIFACTS_SHA256} ${RUNTIME}/ARTIFACTS.json"; do
  printf '%s\n' "${record}" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${SOURCE}" "${RUNTIME}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for path in sys.argv[1:]:
    verify_frozen_tree(Path(path))
PY
module purge
module load python/3.12.3-GCCcore-13.3.0
PYTHONPATH="${SITE_PACKAGES}:${ROOT}" python -m unittest -v \
  tests.unit.test_seurat_wnn_bridge_prepare_smoke \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${SITE_PACKAGES}" python \
  "${ROOT}/scripts/seurat_wnn_bridge_prepare_smoke.py" \
  --inputs "${SOURCE}/prepared/model_inputs" \
  --output "${STAGE}/prepared" \
  > "${STAGE}/validation/prepare.stdout.json"
sha256sum \
  "${ROOT}/scripts/seurat_wnn_bridge_prepare_smoke.py" \
  "${ROOT}/tests/unit/test_seurat_wnn_bridge_prepare_smoke.py" \
  "${ROOT}/slurm/seurat_wnn_bridge_prepare_smoke.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "prepared/receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["pairing_topology"] != "same_nucleus"
    or receipt["outer_unit"] != "donor"
    or set(receipt["folds"]) != set(map(str, range(5)))
    or not receipt["training_pairing_available"]
    or not receipt["query_modality_ids_disjoint"]
    or receipt["hidden_pair_map_read"]
    or receipt["outcomes_read"]
    or receipt["test_outcomes_read"]
    or receipt["rna_to_atac_prediction_executed"]
    or receipt["sealed_rna_conditioned_atac_eligible"]
    or receipt["champion_claim_allowed"]
):
    raise SystemExit("Seurat WNN bridge preparation contract differs")
freeze_tree(root, {
    "artifact_class": "seurat_wnn_bridge_prepared_inputs",
    "dataset_id": "gse296875",
    "source_artifacts_sha256": "c79930c961cbbd44d9d198579f8b9443b6758e28cc2f4458be48d884d62a3354",
    "pairing_topology": "same_nucleus",
    "outer_unit": "donor",
    "hidden_pair_map_read": False,
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
printf 'Seurat WNN bridge prepared inputs: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
