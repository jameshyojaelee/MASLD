#!/usr/bin/env bash
#SBATCH --job-name=spfullcomp
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/spatial_full_composition_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/spatial_full_composition_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON_BIN=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
SCRIPT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"
OUTPUT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/spatial-full-composition-sensitivity-candidate-2026-08-11-r2"

export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
cd "$PROJECT_ROOT"

test ! -e "$OUTPUT_ROOT"

"$PYTHON_BIN" -B -m unittest -v \
  Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/tests/test_full_composition_sensitivity.py \
  Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/tests/test_visium_rerun.py

"$PYTHON_BIN" -B "$SCRIPT_ROOT/31_run_full_composition_sensitivity.py" \
  --project-root "$PROJECT_ROOT" \
  --output-root "$OUTPUT_ROOT"

"$PYTHON_BIN" -B "$SCRIPT_ROOT/32_validate_full_composition_sensitivity.py" \
  --project-root "$PROJECT_ROOT" \
  --output-root "$OUTPUT_ROOT"

python3 scripts/manuscript/validate_resource_scope.py
