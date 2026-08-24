#!/usr/bin/env bash
#SBATCH --job-name=masld-gse267145-qc-rights
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
JOIN="${ROOT}/executions/gse267145-authoritative-join-21064930"
MATRIX="${ROOT}/executions/gse267145-matrix-audit-21064333"
MEASUREMENT="${ROOT}/executions/gse267145-rna-measurement-21065097"
REFERENCE="${ROOT}/executions/gse267145-reference-crosswalk-21065154"
TARGET="${ROOT}/executions/gse267145-qc-rights-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite GSE267145 QC/rights audit: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_gse267145_qc_rights
PYTHONPATH="src" python - "${JOIN}" "${MATRIX}" "${MEASUREMENT}" "${REFERENCE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("input frozen-tree verification passed")
PY
PYTHONPATH="src:." python scripts/audit_gse267145_qc_rights.py \
  --join "${JOIN}" \
  --join-artifacts-sha256 e6c539c5fb29cd357dc126073779948c73c27fd219ff7a4d80b3d20a6cddb580 \
  --matrix "${MATRIX}" \
  --matrix-artifacts-sha256 4e0d17f392530f4e6f0b54346d5d9333e450b7606d9698fbf7253a4f2a0e0417 \
  --measurement "${MEASUREMENT}" \
  --measurement-artifacts-sha256 38531e463837cafd2c8e321411048c7b8b493876c3c964c0f50790cac8806f2f \
  --reference "${REFERENCE}" \
  --reference-artifacts-sha256 da491302ea738a92017978518239bc3884a7635d580e01e3bd28cd4f55b630c2 \
  --output "${STAGE}"
sha256sum \
  scripts/audit_gse267145_qc_rights.py \
  tests/unit/test_audit_gse267145_qc_rights.py \
  slurm/audit_gse267145_qc_rights.sbatch \
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
            "artifact_class": "gse267145_qc_rights_audit",
            "cohort_family_ids": ["gse267145_znf469_human_liver"],
            "biological_unit": "participant",
            "participants": 99,
            "assay_native_qc_frozen": True,
            "internal_nonclinical_research_use": True,
            "data_redistribution_allowed": False,
            "open_weight_release_review_required": True,
            "model_training_activated": False,
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'GSE267145 QC/rights audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
