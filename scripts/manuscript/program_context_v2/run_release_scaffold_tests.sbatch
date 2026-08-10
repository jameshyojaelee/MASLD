#!/usr/bin/env bash
#SBATCH --job-name=manuscript
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/manuscript/program_context_v2/tests/manuscript_%j.out
#SBATCH --error=scripts/manuscript/program_context_v2/tests/manuscript_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
TEST_ROOT="$PROJECT_ROOT/scripts/manuscript/program_context_v2"
SPATIAL_PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python3

export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

cd "$PROJECT_ROOT"
"$SPATIAL_PYTHON" -B "$TEST_ROOT/tests/test_release_scaffold.py"
"$SPATIAL_PYTHON" -B -m unittest -v \
  "$TEST_ROOT/coordinator/tests/test_coordinator_contract.py"
/gpfs/commons/home/jameslee/.local/bin/ruff format --check "$TEST_ROOT"
/gpfs/commons/home/jameslee/.local/bin/ruff check "$TEST_ROOT"
