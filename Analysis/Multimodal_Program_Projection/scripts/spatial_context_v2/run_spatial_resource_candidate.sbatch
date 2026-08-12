#!/usr/bin/env bash
#SBATCH --job-name=spres
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/spatial_resource_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/spatial_resource_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON_BIN=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
R_BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
SCRIPT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"
CANDIDATE_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/spatial-resource-candidate-2026-08-11"

export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
cd "$PROJECT_ROOT"

"$PYTHON_BIN" -B "$SCRIPT_ROOT/20_build_resource_registry.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/20_build_resource_registry.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT" --check-only
"$PYTHON_BIN" -B "$SCRIPT_ROOT/21_build_resource_coverage.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/22_build_hmsma_label_blind.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/23_build_resource_effects.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/24_assemble_resource_candidate.py" \
  --candidate-root "$CANDIDATE_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/25_validate_resource_candidate.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"

export MASLD_FIGURE_OUTPUT_ROOT="$CANDIDATE_ROOT/figures"
"$R_BIN" "$PROJECT_ROOT/scripts/figures/fig4b_protein_triage.R"
"$R_BIN" "$PROJECT_ROOT/scripts/figures/composite_mrna_protein.R"
"$R_BIN" "$PROJECT_ROOT/scripts/figures/fig4d_snatac_accessibility.R"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/26_render_resource_figure4.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/27_validate_resource_figure4.py" \
  --candidate-root "$CANDIDATE_ROOT"
