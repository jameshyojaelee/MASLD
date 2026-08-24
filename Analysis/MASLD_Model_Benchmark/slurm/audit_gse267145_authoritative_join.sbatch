#!/usr/bin/env bash
#SBATCH --job-name=masld-gse267145-join
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
META="${ROOT}/executions/geo-multimodal-source-audit-21064310"
MATRIX="${ROOT}/executions/gse267145-matrix-audit-21064333"
TARGET="${ROOT}/executions/gse267145-authoritative-join-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite GSE267145 join: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_gse267145_authoritative_join
PYTHONPATH="src" python - "${META}" "${MATRIX}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("input frozen-tree verification passed")
PY
PYTHONPATH="src:." python scripts/audit_gse267145_authoritative_join.py \
  --metadata-audit "${META}" \
  --metadata-artifacts-sha256 bac433336705b74dcb2a6ef6b6594af10fef362769d564adfea9788fe546e88e \
  --matrix-audit "${MATRIX}" \
  --matrix-artifacts-sha256 4e0d17f392530f4e6f0b54346d5d9333e450b7606d9698fbf7253a4f2a0e0417 \
  --super-series-url https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267145/soft/GSE267145_family.soft.gz \
  --article-xml-url https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11071482/fullTextXML \
  --output "${STAGE}"
sha256sum \
  scripts/audit_gse267145_authoritative_join.py \
  tests/unit/test_audit_gse267145_authoritative_join.py \
  slurm/audit_gse267145_authoritative_join.sbatch \
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
            "artifact_class": "gse267145_authoritative_participant_join",
            "cohort_family_ids": ["gse267145_znf469_human_liver"],
            "biological_unit": "participant",
            "participants": 99,
            "pairing": "same_sample_different_aliquot",
            "participant_join_authoritative": True,
            "activates_topology_only": True,
            "activates_model_training": False,
            "outcome_visibility": "development_metadata",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'GSE267145 authoritative join: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
