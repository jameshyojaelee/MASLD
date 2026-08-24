#!/usr/bin/env bash
#SBATCH --job-name=masld-atac-transport-outcomes
#SBATCH --partition=io
#SBATCH --cpus-per-task=6
#SBATCH --mem=96G
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
SOURCE_AUDIT="${ROOT}/executions/atac-transport-source-audit-21064435"
SOURCE_SHA=ef43d3c3908c5c81554d750a928adf42f4bec0c7a90f519cdaf062e145e89392
MEMBERSHIP="${ROOT}/executions/atac-transport-membership-21064488"
MEMBERSHIP_SHA=32e1d7b3ba147aea6493a0fcdabe862c1eb887a9a8ec84a1d972b5adbc5b9db6
REGION_ROOT="${ROOT}/executions/chrombpnet-hepatocyte-donor0-genomic0-inputs-21064244/regions"
TARGET="${ROOT}/executions/atac-transport-outcomes-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite ATAC transport outcomes: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python - "${SOURCE_AUDIT}" "${SOURCE_SHA}" "${MEMBERSHIP}" "${MEMBERSHIP_SHA}" <<'PY'
from hashlib import sha256
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for root_text, expected in zip(sys.argv[1::2], sys.argv[2::2], strict=True):
    root = Path(root_text).resolve(strict=True)
    observed = sha256((root / "ARTIFACTS.json").read_bytes()).hexdigest()
    if observed != expected:
        raise SystemExit(f"ARTIFACTS SHA-256 differs: {root}")
    verify_frozen_tree(root)
    print(f"verified frozen input: {root}")
PY
gzip -t /gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/genome/GRCh38.p14.genome.fa.gz
printf '%s  %s\n' \
  9489780d014865df158afc650e1f0ccc204b3a9878e147c1aae2ac982df26215 \
  /gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/genome/GRCh38.p14.genome.fa.gz \
  | sha256sum --check --strict

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
PYTHONPATH="src:." python -m unittest -v tests.unit.test_build_atac_transport_outcomes
PYTHONPATH="src:." python scripts/build_atac_transport_outcomes.py \
  --source-audit "${SOURCE_AUDIT}" \
  --source-artifacts-sha256 "${SOURCE_SHA}" \
  --membership "${MEMBERSHIP}" \
  --membership-artifacts-sha256 "${MEMBERSHIP_SHA}" \
  --regions "${REGION_ROOT}/chrombpnet.regions.bed" \
  --regions-sha256 b1620d5ef4b073b09b9b934d0fe0862678a98ec11e0aadc8fc1c4d684513b3f8 \
  --fold-roles "${REGION_ROOT}/chrombpnet.fold.json" \
  --fold-roles-sha256 c23fb8dd996a6fbc3f15d5f5b61c01a7db2da5d63fbf9fc65a1bf0dd1340ebee \
  --reference-fai /gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa.fai \
  --reference-fai-sha256 bad15b8645d9112f42825660beaff387b2190832ff128a4a8624040a10e0a9b0 \
  --gse281-fragment-builder "${PROJECT}/Analysis/ATAC/Human_External/scripts/02_cellranger_atac_gse281367.sbatch" \
  --gse281-fragment-builder-sha256 2bd9f3630e10f45d9ee9163acfe3c07eb24f3c48a8d87065f9dea0a8165526e6 \
  --gse244-fragment-builder "${PROJECT}/Analysis/ATAC/Human_Multiome/scripts/01_align_scatac.sh" \
  --gse244-fragment-builder-sha256 40eaccc6dea666d99c4667bddcbb1dc22e0e2c38972015dccccbd1b7f94e1726 \
  --gse281-fragments "${PROJECT}/Analysis/ATAC/Human_External/cellranger" \
  --gse244-fragments "${PROJECT}/Analysis/ATAC/Human_Multiome/results/fragments" \
  --profile-bins 20 \
  --workers 4 \
  --output "${STAGE}"
sha256sum \
  scripts/build_atac_transport_outcomes.py \
  tests/unit/test_build_atac_transport_outcomes.py \
  slurm/build_atac_transport_outcomes.sbatch \
  "${PROJECT}/Analysis/ATAC/Human_External/scripts/02_cellranger_atac_gse281367.sbatch" \
  "${PROJECT}/Analysis/ATAC/Human_Multiome/scripts/01_align_scatac.sh" \
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
            "artifact_class": "atac_transport_evaluator_outcomes",
            "cohort_family_ids": ["gse244832", "gse281367"],
            "outcome_visibility": "evaluator_only_development_visible",
            "biological_unit": "donor",
            "fragment_membership_verified": True,
            "same_cell_pairing_claimed": False,
            "activates_development_transport": True,
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'ATAC transport outcomes: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
