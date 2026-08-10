#!/usr/bin/env bash
#SBATCH --job-name=spatial
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/spatial_v2maps_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/spatial_v2maps_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON_BIN=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
SCRIPT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"

export MASLD_PROJECT_ROOT="$PROJECT_ROOT"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
unset FIG4_SPATIAL_SMOKE FIG4_SPATIAL_SMOKE_DATASET FIG4_SPATIAL_N_NULL
unset FIG4_SPATIAL_N_SENSITIVITY_NULL FIG4_SPATIAL_GRAPH_DISTANCE_MULTIPLIER
unset FIG4_SPATIAL_MIN_GRAPH_SPOTS

cd "$PROJECT_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/01_v1_preservation.py" check
"$PYTHON_BIN" -B "$SCRIPT_ROOT/18_build_v2_spatial_map_source.py" \
  --project-root "$PROJECT_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/19_validate_v2_spatial_map_source.py" \
  --project-root "$PROJECT_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/01_v1_preservation.py" check
