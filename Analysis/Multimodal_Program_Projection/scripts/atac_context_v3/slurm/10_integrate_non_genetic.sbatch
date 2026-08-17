#!/bin/bash
#SBATCH --job-name=atacint
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_integrate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_integrate_%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PACKAGE=${ROOT}/Analysis/Multimodal_Program_Projection/scripts/atac_context_v3
CANDIDATE=${ROOT}/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1
MAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
eval "$(${MAMBA} shell hook --shell bash)"
set +u
micromamba activate snapatac2
set -u

python -m py_compile \
  "${PACKAGE}/15_integrate_non_genetic.py" \
  "${ROOT}/Analysis/ATAC/Integration/scripts/35_atac_integration.py"
python "${ROOT}/Analysis/ATAC/Integration/scripts/35_atac_integration.py" \
  --dry-run \
  --outdir "${CANDIDATE}/compatibility/default_before"
python "${PACKAGE}/15_integrate_non_genetic.py"
