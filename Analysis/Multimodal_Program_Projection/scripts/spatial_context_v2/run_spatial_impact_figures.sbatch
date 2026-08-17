#!/usr/bin/env bash
#SBATCH --job-name=spatial_impact
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/spatial_impact_figures_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/spatial_impact_figures_%j.err

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"
OUTPUT_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/spatial-impact-figures-candidate-2026-08-12"
PYTHON_BIN="/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python"

export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

cd "${PROJECT_ROOT}"
test ! -e "${OUTPUT_ROOT}"
"${PYTHON_BIN}" -B "${SCRIPT_ROOT}/49_build_spatial_impact_figures.py" \
  --project-root "${PROJECT_ROOT}" \
  --output-root "${OUTPUT_ROOT}"
"${PYTHON_BIN}" -B "${SCRIPT_ROOT}/50_validate_spatial_impact_figures.py" \
  --project-root "${PROJECT_ROOT}" \
  --output-root "${OUTPUT_ROOT}"
python3 scripts/manuscript/validate_resource_scope.py
