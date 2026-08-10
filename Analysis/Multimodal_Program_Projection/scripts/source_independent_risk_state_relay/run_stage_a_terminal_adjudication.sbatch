#!/bin/bash
#SBATCH --job-name=adjudicate
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay/logs/adjudicate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay/logs/adjudicate_%j.err

set -eo pipefail
set +u
eval "$(micromamba shell hook -s bash)"
micromamba activate /gpfs/commons/home/jameslee/micromamba/envs/rnaseq
set -u
PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
SCRIPT_ROOT=${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay
: "${PLAN45_CANDIDATE_ID:?PLAN45_CANDIDATE_ID must name a new terminal adjudication}"
: "${PLAN45_STAGE_A_EXECUTION_ROOT:?PLAN45_STAGE_A_EXECUTION_ROOT is required}"
: "${PLAN45_STAGE_A_UNBLINDING_ROOT:?PLAN45_STAGE_A_UNBLINDING_ROOT is required}"
: "${PLAN45_STAGE_A_TERMINAL_CONTRACT_ROOT:?PLAN45_STAGE_A_TERMINAL_CONTRACT_ROOT is required}"
: "${PLAN45_STAGE_A_TERMINAL_INPUT_ROOT:?PLAN45_STAGE_A_TERMINAL_INPUT_ROOT is required}"
export MASLD_PROJECT_ROOT=${PROJECT_ROOT}
mkdir -p "${SCRIPT_ROOT}/logs"
cd "${SCRIPT_ROOT}"
python 41_adjudicate_stage_a_and_freeze_stage_b.py
python 42_validate_stage_a_terminal_adjudication.py
