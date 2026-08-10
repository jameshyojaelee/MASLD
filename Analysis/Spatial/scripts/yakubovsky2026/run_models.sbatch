#!/usr/bin/env bash
#SBATCH --job-name=bootstrap
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=72:00:00

set -euo pipefail

BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
PYTHON="/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python"
SCRIPT="${BASE}/Analysis/Spatial/scripts/yakubovsky2026/03_fit_binary_lipid.py"
OUT="${BASE}/Analysis/Spatial/candidates/program-context-v2-candidate-2026-08-07/yakubovsky2026"
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1

"${PYTHON}" "${BASE}/Analysis/Spatial/scripts/yakubovsky2026/06_check_stage.py" \
  --stage scoring --output-dir "${OUT}"

if [[ ! -e "${OUT}/model_stage_seal.tsv" ]]; then
  if [[ -e "${OUT}/donor_program_effects.tsv" ]] || [[ -e "${OUT}/model_execution_manifest.tsv" ]]; then
    printf 'ERROR: partial unsealed model stage exists; refusing resume: %s\n' "${OUT}" >&2
    exit 1
  fi
  "${PYTHON}" "${SCRIPT}" \
    --output-dir "${OUT}" \
    --workers "${SLURM_CPUS_PER_TASK:-16}" \
    --n-bootstrap 4999
else
  "${PYTHON}" "${BASE}/Analysis/Spatial/scripts/yakubovsky2026/06_check_stage.py" \
    --stage model --output-dir "${OUT}"
fi
