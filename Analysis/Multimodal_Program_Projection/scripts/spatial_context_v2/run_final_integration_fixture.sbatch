#!/usr/bin/env bash
#SBATCH --job-name=python
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/python_spint_fixture_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/python_spint_fixture_%j.err

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PYTHON_BIN="/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python"
SCRIPT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"
CANDIDATE_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/spatial_context"

export PYTHONNOUSERSITE=1
cd "$PROJECT_ROOT"

"$PYTHON_BIN" -m unittest discover -s "$SCRIPT_ROOT/tests" -p 'test_final_integration.py' -v
"$PYTHON_BIN" "$SCRIPT_ROOT/13_build_final_integration.py" --mode fixture \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
"$PYTHON_BIN" "$SCRIPT_ROOT/14_validate_final_integration.py" --mode fixture \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
"$PYTHON_BIN" "$SCRIPT_ROOT/14_validate_final_integration.py" --mode fixture --check-only \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
