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

if [[ ! -e "${OUT}/terminal_stage_seal.tsv" ]]; then
  if [[ -e "${OUT}/READY" ]] || [[ -e "${OUT}/validation_report.tsv" ]]; then
    printf 'ERROR: partial unsealed terminal stage exists; refusing resume: %s\n' "${OUT}" >&2
    exit 1
  fi
  "${PYTHON}" "${SCRIPTS}/04_validate_bundle.py" --output-dir "${OUT}"
else
  "${PYTHON}" "${SCRIPTS}/06_check_stage.py" --stage terminal --output-dir "${OUT}"
fi

if [[ ! -e "${OUT}/plan13_adapter/PLAN13_ADAPTER_READY" ]]; then
  if [[ -d "${OUT}/plan13_adapter" ]] && find "${OUT}/plan13_adapter" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
    printf 'ERROR: partial unsealed Plan 13 adapter exists; refusing resume: %s\n' "${OUT}/plan13_adapter" >&2
    exit 1
  fi
  "${PYTHON}" "${SCRIPTS}/05_build_plan13_adapter.py" --write --native-dir "${OUT}"
else
  "${PYTHON}" "${SCRIPTS}/05_build_plan13_adapter.py" --check-only --native-dir "${OUT}"
fi
