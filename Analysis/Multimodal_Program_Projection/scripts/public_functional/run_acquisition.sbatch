#!/bin/bash
#SBATCH --job-name=publicdata
#SBATCH --partition=io
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/logs/acquire_%A_%a.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/logs/acquire_%A_%a.err

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/public_functional"
DATASETS=(GSE200418 GSE207889 GSE253380 GSE106737 CLCC1)

mkdir -p "${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/logs"
cd "${PROJECT_ROOT}"
python3 "${SCRIPT_ROOT}/02_acquire_public_sources.py" --dataset "${DATASETS[${SLURM_ARRAY_TASK_ID}]}"
