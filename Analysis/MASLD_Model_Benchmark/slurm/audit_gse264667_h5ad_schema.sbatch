#!/usr/bin/env bash
#SBATCH --job-name=masld-gse264667-schema
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HDF5_USE_FILE_LOCKING=TRUE
export OMP_NUM_THREADS=1

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT="${PROJECT}/Analysis/MASLD_Model_Benchmark"
INPUT="${ROOT}/executions/gene-perturbation-training-gse264667_hepg2-21067795_0"
SOURCE="${INPUT}/raw/GSE264667_hepg2_raw_singlecell_01.h5ad"
PYTHON="${PROJECT}/.envs/perturbation_gears/bin/python"
TARGET="${ROOT}/executions/gse264667-h5ad-schema-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ ! -x "${PYTHON}" || -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'GSE264667 H5AD schema precondition differs\n' >&2
  exit 1
fi
cd "${ROOT}"
(
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python - "${INPUT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]))
print("source frozen-tree verification passed")
PY
)
"${PYTHON}" scripts/audit_gse264667_h5ad_schema.py \
  --source "${SOURCE}" \
  --output "${STAGE}"
sha256sum \
  scripts/audit_gse264667_h5ad_schema.py \
  slurm/audit_gse264667_h5ad_schema.sbatch \
  config/datasets/gse264667.toml \
  "${INPUT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
(
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
print(freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "gse264667_h5ad_schema_audit",
    "dataset_values_read": False,
    "expression_values_read": False,
    "external_or_soft_links_allowed": False,
    "biological_unit_count": 1,
    "gse313774_accessed": False,
}))
PY
)
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
