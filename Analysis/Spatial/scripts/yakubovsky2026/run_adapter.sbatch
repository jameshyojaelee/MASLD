#!/usr/bin/env bash
#SBATCH --job-name=h5py
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00

set -euo pipefail

BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
PYTHON="/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python"
SCRIPTS="${BASE}/Analysis/Spatial/scripts/yakubovsky2026"
OUT="${BASE}/Analysis/Spatial/candidates/program-context-v2-candidate-2026-08-07/yakubovsky2026"
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "${OUT}/logs"

if [[ ! -e "${OUT}/freeze_stage_seal.tsv" ]]; then
  if compgen -G "${OUT}/analysis_*" >/dev/null || [[ -e "${OUT}/registry_manifest.tsv" ]] || [[ -e "${OUT}/freeze_execution_manifest.tsv" ]]; then
    printf 'ERROR: partial unsealed freeze stage exists; refusing resume: %s\n' "${OUT}" >&2
    exit 1
  fi
  if [[ ! -e "${OUT}/gene_axis_verification.tsv" ]]; then
    "${PYTHON}" "${SCRIPTS}/07_verify_gene_axis.py" \
      --write --output "${OUT}/gene_axis_verification.tsv"
  else
    "${PYTHON}" "${SCRIPTS}/07_verify_gene_axis.py" \
      --check --output "${OUT}/gene_axis_verification.tsv"
  fi
  "${PYTHON}" "${SCRIPTS}/00_freeze_plan11.py" --write --output-dir "${OUT}"
else
  "${PYTHON}" "${SCRIPTS}/06_check_stage.py" --stage freeze --output-dir "${OUT}"
fi

if [[ ! -e "${OUT}/adapter_stage_seal.tsv" ]]; then
  if [[ -e "${OUT}/yakubovsky_human_adapter.h5ad" ]] || [[ -e "${OUT}/adapter_execution_manifest.tsv" ]]; then
    printf 'ERROR: partial unsealed adapter stage exists; refusing resume: %s\n' "${OUT}" >&2
    exit 1
  fi
  "${PYTHON}" "${SCRIPTS}/01_build_adapter.py" --output-dir "${OUT}"
else
  "${PYTHON}" "${SCRIPTS}/06_check_stage.py" --stage adapter --output-dir "${OUT}"
fi

if [[ ! -e "${OUT}/scoring_stage_seal.tsv" ]]; then
  if [[ -e "${OUT}/per_sample_program_scores.tsv" ]] || [[ -e "${OUT}/scoring_execution_manifest.tsv" ]]; then
    printf 'ERROR: partial unsealed scoring stage exists; refusing resume: %s\n' "${OUT}" >&2
    exit 1
  fi
  "${PYTHON}" "${SCRIPTS}/02_score_programs.py" --output-dir "${OUT}"
else
  "${PYTHON}" "${SCRIPTS}/06_check_stage.py" --stage scoring --output-dir "${OUT}"
fi
