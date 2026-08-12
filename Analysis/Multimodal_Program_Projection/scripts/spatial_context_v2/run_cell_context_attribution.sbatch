#!/usr/bin/env bash
#SBATCH --job-name=spcellctx
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/spatial_cell_context_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/spatial_cell_context_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON_BIN=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
SCRIPT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"
OUTPUT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/spatial-cell-context-attribution-candidate-2026-08-11-r2"

export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
cd "$PROJECT_ROOT"

test ! -e "$OUTPUT_ROOT"

"$PYTHON_BIN" -B "$SCRIPT_ROOT/35_build_cell_context_attribution.py" \
  --project-root "$PROJECT_ROOT" \
  --output-root "$OUTPUT_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/36_validate_cell_context_attribution.py" \
  --output-root "$OUTPUT_ROOT"

python3 scripts/manuscript/validate_resource_scope.py
