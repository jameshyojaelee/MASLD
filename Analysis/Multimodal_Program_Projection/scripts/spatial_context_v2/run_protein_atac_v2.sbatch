#!/usr/bin/env bash
#SBATCH --job-name=ATAC
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/ATAC_spint05_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/ATAC_spint05_%j.err

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PYTHON_BIN="/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python"
RSCRIPT_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
SCRIPT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"
CANDIDATE_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/spatial_context"

export MASLD_PROJECT_ROOT="$PROJECT_ROOT"
export SPATIAL_CONTEXT_CANDIDATE_ROOT="$CANDIDATE_ROOT"
export PYTHONNOUSERSITE=1
cd "$PROJECT_ROOT"

"$PYTHON_BIN" -m unittest discover -s "$SCRIPT_ROOT/tests" -p 'test_protein_atac_adapters.py' -v
"$PYTHON_BIN" "$SCRIPT_ROOT/01_v1_preservation.py" check \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
if [[ -d "$CANDIDATE_ROOT/protein_atac_native/atac_v2" ]]; then
  echo "Reusing sealed bounded ATAC native output; downstream validation will rehash every input and producer."
else
  "$RSCRIPT_BIN" "$SCRIPT_ROOT/10_run_atac_v2_projection.R"
fi
"$PYTHON_BIN" "$SCRIPT_ROOT/11_build_protein_atac_adapters.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
"$PYTHON_BIN" "$SCRIPT_ROOT/12_validate_protein_atac_adapters.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
"$PYTHON_BIN" "$SCRIPT_ROOT/12_validate_protein_atac_adapters.py" --check-only \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
"$PYTHON_BIN" "$SCRIPT_ROOT/01_v1_preservation.py" check \
  --project-root "$PROJECT_ROOT" --candidate-root "$CANDIDATE_ROOT"
