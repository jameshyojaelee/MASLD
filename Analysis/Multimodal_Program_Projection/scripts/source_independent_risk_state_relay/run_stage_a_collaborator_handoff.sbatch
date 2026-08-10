#!/bin/bash
#SBATCH --job-name=collabkit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay/logs/collabkit_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay/logs/collabkit_%j.err

set -eo pipefail
set +u
eval "$(micromamba shell hook -s bash)"
micromamba activate /gpfs/commons/home/jameslee/micromamba/envs/rnaseq
set -u

PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
SCRIPT_ROOT=${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay
: "${PLAN45_CANDIDATE_ID:?PLAN45_CANDIDATE_ID must name a new collaborator handoff}"
: "${PLAN45_STAGE_A_TEMPLATE_ROOT:?PLAN45_STAGE_A_TEMPLATE_ROOT is required}"
: "${PLAN45_TARGET_ADJUDICATION_CONTRACT_ROOT:?PLAN45_TARGET_ADJUDICATION_CONTRACT_ROOT is required}"
: "${PLAN45_STAGE_A_EXECUTION_CONTRACT_ROOT:?PLAN45_STAGE_A_EXECUTION_CONTRACT_ROOT is required}"
: "${PLAN45_STAGE_A_QC_CONTRACT_ROOT:?PLAN45_STAGE_A_QC_CONTRACT_ROOT is required}"
: "${PLAN45_STAGE_A_TERMINAL_CONTRACT_ROOT:?PLAN45_STAGE_A_TERMINAL_CONTRACT_ROOT is required}"
export MASLD_PROJECT_ROOT=${PROJECT_ROOT}
mkdir -p "${SCRIPT_ROOT}/logs"
cd "${SCRIPT_ROOT}"
python 43_build_stage_a_collaborator_handoff.py
python 44_validate_stage_a_collaborator_handoff.py
