#!/usr/bin/env bash
#SBATCH --job-name=spatial
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/spatial_semanticv2_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/spatial_semanticv2_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON_BIN=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
SCRIPT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"
RELEASE_ID=program-context-v2-candidate-2026-08-07
SOURCE_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/$RELEASE_ID/spatial_context"
TARGET_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/candidates/$RELEASE_ID/spatial_context_semantic_v2_2026-08-08"
YAK_NATIVE_ROOT="$PROJECT_ROOT/Analysis/Spatial/candidates/$RELEASE_ID/yakubovsky2026"
YAK_ADAPTER_ROOT="$PROJECT_ROOT/Analysis/Spatial/candidates/$RELEASE_ID/semantic_v2_2026-08-08/yakubovsky_plan13_adapter"

export MASLD_PROJECT_ROOT="$PROJECT_ROOT"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

cd "$PROJECT_ROOT"

"$PYTHON_BIN" -B "$SCRIPT_ROOT/16_prepare_semantic_v2.py" \
  --project-root "$PROJECT_ROOT" \
  --source-root "$SOURCE_ROOT" \
  --target-root "$TARGET_ROOT" \
  --yak-native-root "$YAK_NATIVE_ROOT"

"$PYTHON_BIN" -B "$SCRIPT_ROOT/02_build_synthetic_contract.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$TARGET_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/03_validate_contract.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$TARGET_ROOT"

"$PYTHON_BIN" -B "$SCRIPT_ROOT/08_build_real_adapters.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$TARGET_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/09_validate_real_adapters.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$TARGET_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/09_validate_real_adapters.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$TARGET_ROOT" --check-only

"$PYTHON_BIN" -B "$SCRIPT_ROOT/11_build_protein_atac_adapters.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$TARGET_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/12_validate_protein_atac_adapters.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$TARGET_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/12_validate_protein_atac_adapters.py" \
  --project-root "$PROJECT_ROOT" --candidate-root "$TARGET_ROOT" --check-only

"$PYTHON_BIN" -B \
  "$PROJECT_ROOT/Analysis/Spatial/scripts/yakubovsky2026/09_build_plan13_adapter_semantic_v2.py" \
  --base "$PROJECT_ROOT" \
  --native-dir "$YAK_NATIVE_ROOT" \
  --output-dir "$YAK_ADAPTER_ROOT" \
  --write

"$PYTHON_BIN" -B \
  "$PROJECT_ROOT/Analysis/Spatial/scripts/yakubovsky2026/09_build_plan13_adapter_semantic_v2.py" \
  --base "$PROJECT_ROOT" \
  --native-dir "$YAK_NATIVE_ROOT" \
  --output-dir "$YAK_ADAPTER_ROOT" \
  --check-only

"$PYTHON_BIN" -B "$SCRIPT_ROOT/15_freeze_final_integration.py" write \
  --project-root "$PROJECT_ROOT" --candidate-root "$TARGET_ROOT"

"$PYTHON_BIN" -B "$SCRIPT_ROOT/13_build_final_integration.py" --mode fixture \
  --project-root "$PROJECT_ROOT" --candidate-root "$TARGET_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/14_validate_final_integration.py" --mode fixture \
  --project-root "$PROJECT_ROOT" --candidate-root "$TARGET_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/14_validate_final_integration.py" --mode fixture --check-only \
  --project-root "$PROJECT_ROOT" --candidate-root "$TARGET_ROOT"

"$PYTHON_BIN" -B "$SCRIPT_ROOT/13_build_final_integration.py" --mode real \
  --project-root "$PROJECT_ROOT" \
  --candidate-root "$TARGET_ROOT" \
  --yak-native-root "$YAK_NATIVE_ROOT" \
  --yak-adapter-root "$YAK_ADAPTER_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/14_validate_final_integration.py" --mode real \
  --project-root "$PROJECT_ROOT" --candidate-root "$TARGET_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/14_validate_final_integration.py" --mode real --check-only \
  --project-root "$PROJECT_ROOT" --candidate-root "$TARGET_ROOT"

"$PYTHON_BIN" -B "$SCRIPT_ROOT/15_freeze_final_integration.py" check \
  --project-root "$PROJECT_ROOT" --candidate-root "$TARGET_ROOT"
"$PYTHON_BIN" -B "$SCRIPT_ROOT/17_seal_semantic_v2.py" \
  --project-root "$PROJECT_ROOT" \
  --target-root "$TARGET_ROOT" \
  --yak-adapter-root "$YAK_ADAPTER_ROOT"
