#!/usr/bin/env bash
#SBATCH --job-name=masld-gse105127-join-rights
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
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
TARGET="${ROOT}/executions/gse105127-authoritative-join-rights-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite GSE105127 join/rights audit: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_gse105127_join_rights
PYTHONPATH="src" python - "${METADATA}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
print("metadata frozen-tree verification passed")
PY
PYTHONPATH="src:." python scripts/audit_gse105127_join_rights.py \
  --metadata "${METADATA}" \
  --metadata-artifacts-sha256 bac433336705b74dcb2a6ef6b6594af10fef362769d564adfea9788fe546e88e \
  --output "${STAGE}"
sha256sum \
  scripts/audit_gse105127_join_rights.py \
  tests/unit/test_audit_gse105127_join_rights.py \
  slurm/audit_gse105127_join_rights.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "gse105127_authoritative_join_rights_audit",
            "cohort_family_ids": ["gse105127_zonated_rna_rrbs"],
            "biological_unit": "participant",
            "participants": 19,
            "assay_records": 114,
            "pairing": "adjacent_section",
            "participant_safe_outer_folds": True,
            "internal_nonclinical_research_use": True,
            "data_redistribution_allowed": False,
            "model_training_activated": False,
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'GSE105127 authoritative join/rights audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
