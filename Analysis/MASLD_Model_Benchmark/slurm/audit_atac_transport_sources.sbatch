#!/usr/bin/env bash
#SBATCH --job-name=masld-atac-transport-audit
#SBATCH --partition=io
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT="${PROJECT}/Analysis/MASLD_Model_Benchmark"
TARGET="${ROOT}/executions/atac-transport-source-audit-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite ATAC transport audit: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_atac_transport_sources

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
PYTHONPATH="src:." python scripts/audit_atac_transport_sources.py \
  --gse281-h5ad "${PROJECT}/Analysis/ATAC/Human_External/snapatac2_fast/snapatac2_processed.h5ad" \
  --gse281-meta "${PROJECT}/data/GSE281367/metadata/donor_pairing.csv" \
  --gse281-fragments "${PROJECT}/Analysis/ATAC/Human_External/cellranger" \
  --gse281-pseudobulk "${PROJECT}/Analysis/ATAC/Human_External/pseudobulk" \
  --gse244-h5ad "${PROJECT}/Analysis/ATAC/Human_Multiome/results/label_transfer/snapatac2_relabeled.h5ad" \
  --gse244-meta "${PROJECT}/Analysis/ATAC/Human_Multiome/metadata/donor_metadata_curated.tsv" \
  --gse244-fragments "${PROJECT}/Analysis/ATAC/Human_Multiome/results/fragments" \
  --gse244-pseudobulk "${PROJECT}/Analysis/ATAC/Human_Multiome/results/snapatac2" \
  --hash-workers 8 \
  --output "${STAGE}"
sha256sum \
  scripts/audit_atac_transport_sources.py \
  tests/unit/test_audit_atac_transport_sources.py \
  slurm/audit_atac_transport_sources.sbatch \
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
            "artifact_class": "atac_transport_source_admission_audit",
            "cohort_family_ids": ["gse244832", "gse281367"],
            "outcome_visibility": "development_visible",
            "activates_dataset": False,
            "legacy_pseudobulk_is_benchmark_authority": False,
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'ATAC transport source audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
