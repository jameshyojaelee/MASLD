#!/usr/bin/env bash
#SBATCH --job-name=masld-atac-transport-membership
#SBATCH --partition=io
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT="${PROJECT}/Analysis/MASLD_Model_Benchmark"
SOURCE_AUDIT="${ROOT}/executions/atac-transport-source-audit-21064435"
SOURCE_ARTIFACTS_SHA256=ef43d3c3908c5c81554d750a928adf42f4bec0c7a90f519cdaf062e145e89392
TARGET="${ROOT}/executions/atac-transport-membership-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite ATAC transport membership: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_atac_transport_sources \
  tests.unit.test_build_atac_transport_membership
PYTHONPATH="src:." python - "${SOURCE_AUDIT}" "${SOURCE_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import hashlib
import sys

from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
if hashlib.sha256((root / "ARTIFACTS.json").read_bytes()).hexdigest() != sys.argv[2]:
    raise SystemExit("source-audit ARTIFACTS SHA-256 differs")
verify_frozen_tree(root)
print("source-audit frozen tree verified")
PY

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
PYTHONPATH="src:." python scripts/build_atac_transport_membership.py \
  --source-audit "${SOURCE_AUDIT}" \
  --source-artifacts-sha256 "${SOURCE_ARTIFACTS_SHA256}" \
  --gse281-merged "${PROJECT}/Analysis/ATAC/Human_External/snapatac2_fast/snapatac2_processed.h5ad" \
  --gse281-per-donor "${PROJECT}/Analysis/ATAC/Human_External/snapatac2_fast/per_donor" \
  --gse281-meta "${PROJECT}/data/GSE281367/metadata/donor_pairing.csv" \
  --gse244-merged "${PROJECT}/Analysis/ATAC/Human_Multiome/results/label_transfer/snapatac2_relabeled.h5ad" \
  --gse244-per-donor "${PROJECT}/Analysis/ATAC/Human_Multiome/results/snapatac2/per_donor" \
  --gse244-meta "${PROJECT}/Analysis/ATAC/Human_Multiome/metadata/donor_metadata_curated.tsv" \
  --hash-workers 8 \
  --output "${STAGE}"
sha256sum \
  scripts/audit_atac_transport_sources.py \
  scripts/build_atac_transport_membership.py \
  tests/unit/test_audit_atac_transport_sources.py \
  tests/unit/test_build_atac_transport_membership.py \
  slurm/build_atac_transport_membership.sbatch \
  > "${STAGE}/source.sha256"
micromamba deactivate
PYTHONPATH="src:." python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "atac_transport_cell_membership",
            "cohort_family_ids": ["gse244832", "gse281367"],
            "source_audit_artifacts_sha256": "ef43d3c3908c5c81554d750a928adf42f4bec0c7a90f519cdaf062e145e89392",
            "outcome_visibility": "development_visible",
            "fragment_membership_verified": False,
            "activates_dataset": False,
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'ATAC transport membership: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
