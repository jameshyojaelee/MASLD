#!/usr/bin/env bash
#SBATCH --job-name=masld-geo-multimodal-audit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/geo-multimodal-source-audit-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite source audit: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_geo_multimodal_sources \
  tests.contract.test_cross_cohort_plan
PYTHONPATH="src:." python scripts/audit_geo_multimodal_sources.py \
  --output "${STAGE}"
sha256sum \
  scripts/audit_geo_multimodal_sources.py \
  tests/unit/test_audit_geo_multimodal_sources.py \
  tests/contract/test_cross_cohort_plan.py \
  config/datasets/gse267145_znf469_human_liver.toml \
  config/datasets/gse105127_zonated_rna_rrbs.toml \
  config/evaluation/cross_cohort_expansion.toml \
  config/evaluation/family_native_tournament.toml \
  slurm/audit_geo_multimodal_sources.sbatch \
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
            "artifact_class": "geo_multimodal_source_admission_audit",
            "cohort_family_ids": [
                "gse267145_znf469_human_liver",
                "gse105127_zonated_rna_rrbs",
            ],
            "outcome_visibility": "development_metadata",
            "activates_dataset": False,
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'GEO multimodal source audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
