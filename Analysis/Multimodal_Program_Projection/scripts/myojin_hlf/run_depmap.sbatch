#!/usr/bin/env bash
#SBATCH --job-name=depmap
#SBATCH --partition=io
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf/logs/depmap_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf/logs/depmap_%j.err

set -euo pipefail
PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
bash "${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/myojin_hlf/acquire_depmap_24q4.sh"
