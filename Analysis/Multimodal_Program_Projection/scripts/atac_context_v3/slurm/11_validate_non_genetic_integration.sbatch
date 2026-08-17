#!/bin/bash
#SBATCH --job-name=atacintv
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_integration_validate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_integration_validate_%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PACKAGE=${ROOT}/Analysis/Multimodal_Program_Projection/scripts/atac_context_v3
CANDIDATE=${ROOT}/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1
INTEGRATOR=${ROOT}/Analysis/ATAC/Integration/scripts/35_atac_integration.py
ATLAS=${ROOT}/RNA-seq/results/multi_evidence/multi_evidence_atlas.csv
MAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
eval "$(${MAMBA} shell hook --shell bash)"
set +u
micromamba activate snapatac2
set -u

canonical_before=$(sha256sum "${ATLAS}" | awk '{print $1}')
python -m py_compile "${INTEGRATOR}" "${PACKAGE}/16_validate_non_genetic_integration.py"
python "${INTEGRATOR}" \
  --dry-run \
  --outdir "${CANDIDATE}/compatibility/default_after"
python "${INTEGRATOR}" \
  --dry-run \
  --atac-v3-root "${CANDIDATE}" \
  --outdir "${CANDIDATE}/compatibility/with_v3"
if python "${INTEGRATOR}" \
  --atac-v3-root "${CANDIDATE}" \
  --outdir "${CANDIDATE}/compatibility/unsafe_write_probe"; then
  echo "ERROR: candidate mode allowed a non-dry-run write" >&2
  exit 1
fi
python "${PACKAGE}/16_validate_non_genetic_integration.py" \
  --canonical-before-sha256 "${canonical_before}"
