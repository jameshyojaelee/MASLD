#!/usr/bin/env bash
#SBATCH --job-name=spunitnull
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/spatial_unit_null_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/spatial_unit_null_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON_BIN=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
SCRIPT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"
UPSTREAM_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/spatial_context/native_spatial/v2_candidate"
OUTPUT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/spatial-unit-null-diagnostics-candidate-2026-08-11-r2"

export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
cd "$PROJECT_ROOT"

test -f "$UPSTREAM_ROOT/READY"
test ! -e "$OUTPUT_ROOT"

"$PYTHON_BIN" -B "$SCRIPT_ROOT/33_render_unit_null_diagnostics.py" \
  --upstream-root "$UPSTREAM_ROOT" \
  --output-root "$OUTPUT_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/34_validate_unit_null_diagnostics.py" \
  --upstream-root "$UPSTREAM_ROOT" \
  --output-root "$OUTPUT_ROOT"

python3 scripts/manuscript/validate_resource_scope.py
