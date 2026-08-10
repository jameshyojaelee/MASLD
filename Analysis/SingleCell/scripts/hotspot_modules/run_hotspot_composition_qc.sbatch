#!/usr/bin/env bash
#SBATCH --job-name=limma
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/logs/limma_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/logs/limma_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPT_ROOT="$PROJECT_ROOT/Analysis/SingleCell/scripts/hotspot_modules"
OUTPUT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
RSCRIPT=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript

cd "$PROJECT_ROOT"
"$RSCRIPT" "$SCRIPT_ROOT/517_refit_qc_composition.R" \
  --project-root "$PROJECT_ROOT" \
  --output-root "$OUTPUT_ROOT" \
  --producer-script "$SCRIPT_ROOT/517_refit_qc_composition.R"
"$RSCRIPT" "$SCRIPT_ROOT/518_validate_qc_composition.R" \
  --project-root "$PROJECT_ROOT" \
  --output-root "$OUTPUT_ROOT" \
  --producer-script "$SCRIPT_ROOT/518_validate_qc_composition.R"
