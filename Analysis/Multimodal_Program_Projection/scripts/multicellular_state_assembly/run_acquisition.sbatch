#!/usr/bin/env bash
#SBATCH --job-name=p44acq
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/multicellular-assembly-response-2026-08-09/logs/p44acq_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/multicellular-assembly-response-2026-08-09/logs/p44acq_%j.err

set -euo pipefail
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
mkdir -p "${ROOT}/Analysis/Multimodal_Program_Projection/candidates/multicellular-assembly-response-2026-08-09/logs"
bash "${ROOT}/Analysis/Multimodal_Program_Projection/scripts/multicellular_state_assembly/04_acquire_expression_sources.sh"
