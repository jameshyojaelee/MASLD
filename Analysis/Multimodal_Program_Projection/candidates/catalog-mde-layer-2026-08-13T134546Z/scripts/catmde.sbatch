#!/bin/bash
#SBATCH --job-name=catmde
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/catalog-mde-layer-2026-08-13T134546Z/logs/catmde_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/catalog-mde-layer-2026-08-13T134546Z/logs/catmde_%j.err

set -euo pipefail
export PYTHONNOUSERSITE=1
export MDE_OUT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/catalog-mde-layer-2026-08-13T134546Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python

"$PY" -m pip freeze > "$MDE_OUT/env/pip_freeze.txt"
"$PY" -c "import sys; print(sys.version)" > "$MDE_OUT/env/python_version.txt"
"$PY" "$MDE_OUT/scripts/01_build_mde_layer.py"
echo "DONE" > "$MDE_OUT/BUILD_COMPLETE"
