#!/bin/bash
# Wrapper: runs the sealed 09_genetics_validate.sbatch unmodified, with memory
# raised from its own 32G header to 128G after job 20918763 was killed
# OUT_OF_MEMORY at 30.7GB, and on --qos=interactive per user instruction.
#SBATCH --job-name=atacgenval
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_genvalidate_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_genvalidate_%j.log

set -euo pipefail
SEALED=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/slurm/09_genetics_validate.sbatch
EXPECT=f44e1ba81518e5fb6fd3c8d4e552a5b51f96b3300beaeebf31b981296633f62d
ACTUAL=$(sha256sum "${SEALED}" | cut -d' ' -f1)
if [[ "${ACTUAL}" != "${EXPECT}" ]]; then
  echo "sealed genetics-validate script hash mismatch: ${ACTUAL}" >&2
  exit 1
fi
exec bash "${SEALED}"
