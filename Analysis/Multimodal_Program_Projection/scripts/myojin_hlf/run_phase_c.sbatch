#!/usr/bin/env bash
#SBATCH --job-name=permutation
#SBATCH --partition=cpu
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf/logs/permutation_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf/logs/permutation_%j.err

set -euo pipefail
PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
CANDIDATE_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf"
SCRIPT_DIR="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/myojin_hlf"

export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

micromamba run -n rnaseq python -B "${SCRIPT_DIR}/10_execute_one_pass.py" --outdir "${CANDIDATE_ROOT}"
micromamba run -n rnaseq python -B "${SCRIPT_DIR}/11_validate_phase_c.py" --outdir "${CANDIDATE_ROOT}"
