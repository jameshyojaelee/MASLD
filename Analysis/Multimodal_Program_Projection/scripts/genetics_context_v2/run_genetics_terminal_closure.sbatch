#!/bin/bash
#SBATCH --job-name=eqtl
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/genetics_context/logs/eqtl_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/genetics_context/logs/eqtl_%j.err

set -euo pipefail

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
SCRIPT_DIR="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/scripts/genetics_context_v2"
CANDIDATE_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/genetics_context"
PYTHON_BIN="${PYTHON_BIN:-/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python3}"

test -d "$CANDIDATE_ROOT/logs" || {
  echo "Missing pre-created SLURM log directory: $CANDIDATE_ROOT/logs" >&2
  exit 2
}
test -x "$PYTHON_BIN" || {
  echo "Missing executable Python: $PYTHON_BIN" >&2
  exit 2
}

export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1

"$PYTHON_BIN" "$SCRIPT_DIR/09_close_genetics_terminal_gate.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
"$PYTHON_BIN" "$SCRIPT_DIR/10_validate_genetics_terminal_closure.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"

test -s "$CANDIDATE_ROOT/GEN_TERMINAL_CLOSURE_READY"
echo "GEN terminal closure validated; historical preflight gate was not overwritten"
