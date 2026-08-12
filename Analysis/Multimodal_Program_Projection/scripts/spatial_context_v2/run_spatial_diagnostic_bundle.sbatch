#!/usr/bin/env bash
#SBATCH --job-name=spdiagbundle
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/spatial_diagnostic_bundle_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/spatial_diagnostic_bundle_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON_BIN=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
SCRIPT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"
OUTPUT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/spatial-diagnostic-bundle-candidate-2026-08-11"

export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
cd "$PROJECT_ROOT"

test ! -e "$OUTPUT_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/45_build_spatial_diagnostic_bundle.py" \
  --project-root "$PROJECT_ROOT" --output-root "$OUTPUT_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/46_validate_spatial_diagnostic_bundle.py" \
  --project-root "$PROJECT_ROOT" --output-root "$OUTPUT_ROOT"
python3 scripts/manuscript/validate_resource_scope.py
