#!/bin/bash
#SBATCH --job-name=pairorient
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay/logs/pairorient_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay/logs/pairorient_%A_%a.err

set -eo pipefail
set +u
eval "$(micromamba shell hook -s bash)"
micromamba activate /gpfs/commons/home/jameslee/micromamba/envs/rnaseq
set -u

PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
SCRIPT_ROOT=${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay
mkdir -p "${SCRIPT_ROOT}/logs"

: "${PLAN45_CANDIDATE_ID:?export the sealed orientation candidate ID}"
: "${SLURM_ARRAY_TASK_ID:?submit as a zero-based array over every orientation-worklist row}"
export MASLD_PROJECT_ROOT=${PROJECT_ROOT}
export PLAN45_ORIENTATION_INDEX=${SLURM_ARRAY_TASK_ID}
export LD_PANEL=polyfun
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

Rscript "${SCRIPT_ROOT}/06_targeted_susie_pair_export.R"
