#!/bin/bash
#SBATCH --job-name=spatial
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/spatial_v1semantic_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/spatial_v1semantic_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
VALIDATOR="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/06b_validate_visium_v1_semantic_equivalence.py"

export MASLD_PROJECT_ROOT="$PROJECT_ROOT"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

cd "$PROJECT_ROOT"
micromamba run -n spatial python -B "$VALIDATOR" \
  --check-only --project-root "$PROJECT_ROOT"
micromamba run -n spatial python -B "$VALIDATOR" \
  --project-root "$PROJECT_ROOT"
