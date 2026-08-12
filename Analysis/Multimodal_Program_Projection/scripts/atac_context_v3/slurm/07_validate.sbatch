#!/bin/bash
#SBATCH --job-name=atacvalid
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_validate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_validate_%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PACKAGE=${ROOT}/Analysis/Multimodal_Program_Projection/scripts/atac_context_v3
MAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
eval "$(${MAMBA} shell hook --shell bash)"
micromamba activate snapatac2
python "${PACKAGE}/10a_validate_exact_recount.py"
python "${PACKAGE}/10_validate_non_genetic.py"
python "${PACKAGE}/07_prepare_genetic_replay.py"
