#!/bin/bash
#SBATCH --job-name=pubfreeze
#SBATCH --partition=cpu
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/logs/freeze_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/logs/freeze_%j.err

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/public_functional"
mkdir -p "${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/logs"
cd "${PROJECT_ROOT}"
python3 "${SCRIPT_ROOT}/01_validate_seal.py"
