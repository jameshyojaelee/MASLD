#!/usr/bin/env bash
#SBATCH --job-name=masld-gse105127-rrbs-qc-r2
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
JOIN="${ROOT}/executions/gse105127-authoritative-join-rights-21065552"
REUSE="${ROOT}/executions/gse105127-rrbs-source-qc-21065724.failed/raw"
TARGET="${ROOT}/executions/gse105127-rrbs-source-qc-r2-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed GSE105127 RRBS r2 stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE105127 RRBS r2 audit: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_acquire_audit_gse105127_rrbs
PYTHONPATH="src" python - "${JOIN}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
print("authoritative join frozen-tree verification passed")
PY
test "$(find "${REUSE}" -maxdepth 1 -type f -name '*.bed.gz' | wc -l)" -eq 57
PYTHONPATH="src:." python scripts/acquire_audit_gse105127_rrbs.py \
  --join "${JOIN}" \
  --join-artifacts-sha256 574eeee1119054a2d7b55066804accf29fb45f23b591dd26f9d61d907b27a75f \
  --reuse-raw "${REUSE}" \
  --output "${STAGE}"
sha256sum \
  scripts/acquire_audit_gse105127_rrbs.py \
  tests/unit/test_acquire_audit_gse105127_rrbs.py \
  slurm/acquire_audit_gse105127_rrbs_reuse.sbatch \
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
            "artifact_class": "gse105127_rrbs_source_qc",
            "cohort_family_ids": ["gse105127_zonated_rna_rrbs"],
            "biological_unit": "participant",
            "participants": 19,
            "rrbs_files": 57,
            "pairing": "adjacent_section",
            "coverage_aware_qc": True,
            "bed_detail_auxiliary_semantics_resolved": False,
            "coordinate_remap_complete": False,
            "model_training_activated": False,
        },
    )
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE105127 RRBS source/QC r2 audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
