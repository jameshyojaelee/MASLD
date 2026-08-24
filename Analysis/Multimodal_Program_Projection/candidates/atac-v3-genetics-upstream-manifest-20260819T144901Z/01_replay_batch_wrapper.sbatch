#!/bin/bash
# Wrapper: creates the node-local /scratch/$USER parent that the sealed
# 08_genetics_batch.sbatch assumes, then runs that sealed script unmodified.
# Resources below are copied verbatim from the sealed script's own header.
#SBATCH --job-name=atacreplay
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_replay_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_replay_%j.err

set -euo pipefail
if [[ $# -ne 1 || ! $1 =~ ^[1-5]$ ]]; then
  echo "usage: sbatch 01_replay_batch_wrapper.sbatch 1|2|3|4|5" >&2
  exit 2
fi
mkdir -p "/scratch/${USER}"
SEALED=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/slurm/08_genetics_batch.sbatch
EXPECT=f1a8d7f01c0058098542d442e7c7e61defe307de751ffafcccf4753bb33e372d
ACTUAL=$(sha256sum "${SEALED}" | cut -d' ' -f1)
if [[ "${ACTUAL}" != "${EXPECT}" ]]; then
  echo "sealed replay script hash mismatch: ${ACTUAL}" >&2
  exit 1
fi
exec bash "${SEALED}" "$1"
