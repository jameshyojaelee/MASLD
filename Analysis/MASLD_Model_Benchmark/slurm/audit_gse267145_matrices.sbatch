#!/usr/bin/env bash
#SBATCH --job-name=masld-gse267145-matrix-audit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
METADATA="${ROOT}/executions/geo-multimodal-source-audit-21064310"
METADATA_ARTIFACTS_SHA256=bac433336705b74dcb2a6ef6b6594af10fef362769d564adfea9788fe546e88e
TARGET="${ROOT}/executions/gse267145-matrix-audit-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite matrix audit: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_gse267145_matrices
PYTHONPATH="src:." python scripts/audit_gse267145_matrices.py \
  --metadata-audit "${METADATA}" \
  --metadata-artifacts-sha256 "${METADATA_ARTIFACTS_SHA256}" \
  --output "${STAGE}"
sha256sum \
  scripts/audit_gse267145_matrices.py \
  tests/unit/test_audit_gse267145_matrices.py \
  slurm/audit_gse267145_matrices.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="src:." python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "gse267145_matrix_admission_audit",
            "cohort_family_id": "gse267145_znf469_human_liver",
            "metadata_artifacts_sha256": "bac433336705b74dcb2a6ef6b6594af10fef362769d564adfea9788fe546e88e",
            "outcome_visibility": "development_matrix",
            "activates_dataset": False,
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'GSE267145 matrix audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
