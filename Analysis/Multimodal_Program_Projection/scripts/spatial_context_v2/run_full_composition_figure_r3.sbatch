#!/usr/bin/env bash
#SBATCH --job-name=spcompf3
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/spatial_full_composition_figure_r3_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/spatial_full_composition_figure_r3_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON_BIN=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
SCRIPT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"
INPUT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/spatial-full-composition-sensitivity-candidate-2026-08-11-r2"
OUTPUT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/spatial-full-composition-figure-candidate-2026-08-11-r3"

export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
cd "$PROJECT_ROOT"

test -f "$INPUT_ROOT/READY"
test ! -e "$OUTPUT_ROOT"

"$PYTHON_BIN" -B "$SCRIPT_ROOT/41_render_full_composition_figure_r3.py" \
  --project-root "$PROJECT_ROOT" \
  --input-root "$INPUT_ROOT" \
  --output-root "$OUTPUT_ROOT"

"$PYTHON_BIN" -B "$SCRIPT_ROOT/42_validate_full_composition_figure_r3.py" \
  --project-root "$PROJECT_ROOT" \
  --input-root "$INPUT_ROOT" \
  --output-root "$OUTPUT_ROOT"
