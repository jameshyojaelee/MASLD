#!/usr/bin/env bash
#SBATCH --job-name=masld-bulk-source-audit
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/bulk-expansion-source-audit-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite bulk source audit: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_bulk_expansion_sources
PYTHONPATH="src:." python scripts/audit_bulk_expansion_sources.py \
  --output "${STAGE}"
sha256sum \
  scripts/audit_bulk_expansion_sources.py \
  tests/unit/test_audit_bulk_expansion_sources.py \
  slurm/audit_bulk_expansion_sources.sbatch \
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
            "artifact_class": "bulk_expansion_source_audit",
            "cohort_family_ids": [
                "gse268273_imid_masld",
                "gse260666_bulk_rna",
                "gse274114_mash_hbv",
            ],
            "biological_unit": "candidate_participant",
            "official_geo_sample_records": 164,
            "participant_join_authoritative": False,
            "model_training_activated": False,
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'Bulk expansion source audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
