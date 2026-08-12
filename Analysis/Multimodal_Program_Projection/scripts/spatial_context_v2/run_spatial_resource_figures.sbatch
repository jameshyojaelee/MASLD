#!/usr/bin/env bash
#SBATCH --job-name=spfig
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/spatial_resource_figures_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/spatial_resource_figures_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON_BIN=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
R_BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
SCRIPT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"
CANDIDATE_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/spatial-resource-candidate-2026-08-11"

export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export MASLD_PROJECT_ROOT="$PROJECT_ROOT"
export MASLD_FIGURE_OUTPUT_ROOT="$CANDIDATE_ROOT/figures"
cd "$PROJECT_ROOT"

test -f "$CANDIDATE_ROOT/validation/READY"
test ! -e "$CANDIDATE_ROOT/figures/READY"

"$R_BIN" "$PROJECT_ROOT/scripts/figures/fig4b_protein_triage.R"
"$R_BIN" "$PROJECT_ROOT/scripts/figures/composite_mrna_protein.R"
"$R_BIN" "$PROJECT_ROOT/scripts/figures/fig4d_snatac_accessibility.R"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/26_render_resource_figure4.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/27_validate_resource_figure4.py" \
  --candidate-root "$CANDIDATE_ROOT"
