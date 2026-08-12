#!/usr/bin/env bash
#SBATCH --job-name=spverify
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/spatial_final_verification_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/spatial_final_verification_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON_BIN=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
SCRIPT_ROOT="$PROJECT_ROOT/Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"
VERIFY_PYCACHE=$(mktemp -d /tmp/spatial-final-verification-pycache.XXXXXX)
export PYTHONPYCACHEPREFIX="$VERIFY_PYCACHE"
trap 'rm -rf -- "$VERIFY_PYCACHE"' EXIT
cd "$PROJECT_ROOT"

mapfile -t SPATIAL_SCRIPTS < <(find "$SCRIPT_ROOT" -maxdepth 1 -type f -name '*.py' | sort)
"$PYTHON_BIN" -B -m py_compile "${SPATIAL_SCRIPTS[@]}"

"$PYTHON_BIN" -B -m unittest -v \
  Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/tests/test_final_integration.py \
  Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/tests/test_full_composition_sensitivity.py \
  Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/tests/test_protein_atac_adapters.py \
  Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/tests/test_real_adapters.py \
  Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/tests/test_spatial_resource.py \
  Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/tests/test_visium_rerun.py

python3 scripts/manuscript/validate_resource_scope.py

set +e
"$PYTHON_BIN" -B scripts/manuscript/check_release_consistency.py
release_consistency_exit=$?
set -e
echo "RELEASE_CONSISTENCY_EXIT=$release_consistency_exit"
