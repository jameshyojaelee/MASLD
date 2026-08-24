#!/usr/bin/env bash
#SBATCH --job-name=masld-bulk-matrix-audit
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE="${ROOT}/executions/bulk-expansion-source-audit-21065762"
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python
TARGET="${ROOT}/executions/bulk-expansion-matrix-audit-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite bulk matrix audit: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
"${BASE}" -m unittest -v tests.unit.test_acquire_audit_bulk_expansion_matrices
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src" python - "${SOURCE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
print("bulk source frozen-tree verification passed")
PY
PYTHONPATH="src:." "${BASE}" scripts/acquire_audit_bulk_expansion_matrices.py \
  --source "${SOURCE}" \
  --source-artifacts-sha256 883f870e2a76f78c39d257c24f8e795d2e74f1bd32aba785251b52aa7cb8229e \
  --output "${STAGE}"
sha256sum \
  scripts/acquire_audit_bulk_expansion_matrices.py \
  tests/unit/test_acquire_audit_bulk_expansion_matrices.py \
  slurm/acquire_audit_bulk_expansion_matrices.sbatch \
  > "${STAGE}/source.sha256"
"${BASE}" -m pip freeze > "${STAGE}/effective-pip-freeze.txt"
PYTHONPATH="src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "bulk_expansion_matrix_audit",
            "cohort_family_ids": [
                "gse268273_imid_masld",
                "gse260666_bulk_rna",
                "gse274114_mash_hbv",
            ],
            "biological_unit": "candidate_participant",
            "processed_sources": 42,
            "participant_join_authoritative": False,
            "model_training_activated": False,
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'Bulk expansion matrix audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
