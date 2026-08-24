#!/bin/bash
#SBATCH --job-name=atacprep
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/atac-v3-genetics-upstream-manifest-20260819T144901Z/prepare_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/atac-v3-genetics-upstream-manifest-20260819T144901Z/prepare_%j.log

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PACKAGE=${ROOT}/Analysis/Multimodal_Program_Projection/scripts/atac_context_v3
WORK=${ROOT}/Analysis/Multimodal_Program_Projection/candidates/atac-v3-genetics-upstream-manifest-20260819T144901Z
MAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
eval "$(${MAMBA} shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u
cd "${PACKAGE}"
python 07_prepare_genetic_replay.py --promoted-manifest "${WORK}/promoted_coloc_manifest.tsv"
