#!/usr/bin/env bash
#SBATCH --job-name=spdownfig
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/spatial_downstream_figures_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/spatial_downstream_figures_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON_BIN=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
SCRIPT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"
UPSTREAM_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/spatial-resource-candidate-2026-08-11"
OUTPUT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/spatial-figure-integration-candidate-2026-08-11-r2"

export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
cd "$PROJECT_ROOT"

test -f "$UPSTREAM_ROOT/validation/READY"
test -f "$UPSTREAM_ROOT/figures/READY"
test ! -e "$OUTPUT_ROOT"

"$PYTHON_BIN" -B "$SCRIPT_ROOT/29_render_resource_downstream_figures.py" \
  --project-root "$PROJECT_ROOT" \
  --upstream-root "$UPSTREAM_ROOT" \
  --output-root "$OUTPUT_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/30_validate_resource_downstream_figures.py" \
  --upstream-root "$UPSTREAM_ROOT" \
  --output-root "$OUTPUT_ROOT"

python3 scripts/manuscript/validate_resource_scope.py
