#!/usr/bin/env bash
#SBATCH --job-name=python
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00

set -euo pipefail

BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
PYTHON="/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python"
SCRIPTS="${BASE}/Analysis/Spatial/scripts/yakubovsky2026"
OUT="${BASE}/Analysis/Spatial/candidates/program-context-v2-candidate-2026-08-07/yakubovsky2026"
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"

"${PYTHON}" "${SCRIPTS}/06_check_stage.py" --stage model --output-dir "${OUT}"
"${PYTHON}" "${SCRIPTS}/08_validate_bundle_hotfix.py" --output-dir "${OUT}"
"${PYTHON}" "${SCRIPTS}/05_build_plan13_adapter.py" --write --native-dir "${OUT}"
