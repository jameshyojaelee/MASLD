#!/bin/bash
#SBATCH --job-name=coloc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/genetics_context/logs/coloc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/genetics_context/logs/coloc_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPT_DIR="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/scripts/genetics_context_v2"
CANDIDATE_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/genetics_context"

mkdir -p "$CANDIDATE_ROOT/logs"
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1

micromamba run -n rnaseq python3 "$SCRIPT_DIR/01_freeze_and_rederive.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
micromamba run -n rnaseq python3 "$SCRIPT_DIR/02_build_phenotype_registry.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
micromamba run -n rnaseq python3 "$SCRIPT_DIR/03_build_eqtl_observability.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
micromamba run -n rnaseq Rscript "$SCRIPT_DIR/04_power_interface.R" \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
micromamba run -n rnaseq python3 "$SCRIPT_DIR/04_seal_power_stage.py"
micromamba run -n rnaseq python3 "$SCRIPT_DIR/05_validate_preflight.py"

