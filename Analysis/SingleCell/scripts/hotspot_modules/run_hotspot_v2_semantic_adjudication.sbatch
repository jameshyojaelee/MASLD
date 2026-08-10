#!/bin/bash
#SBATCH --job-name=hotspot
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/logs/hotspot_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/logs/hotspot_%j.err

set -euo pipefail

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
SCRIPT_DIR="$PROJECT_ROOT/Analysis/SingleCell/scripts/hotspot_modules"
CANDIDATE_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
PYTHON_BIN="${PYTHON_BIN:-python3}"

test -d "$CANDIDATE_ROOT/logs" || {
  echo "Missing pre-created SLURM log directory: $CANDIDATE_ROOT/logs" >&2
  exit 2
}

export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1

"$PYTHON_BIN" "$SCRIPT_DIR/515_adjudicate_hotspot_v2_semantics.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
"$PYTHON_BIN" "$SCRIPT_DIR/516_validate_hotspot_v2_semantics.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"

test -s "$CANDIDATE_ROOT/SEMANTIC_ADJUDICATION_READY"
echo "Hotspot v2 semantic adjudication validated; frozen registry was not rewritten"
