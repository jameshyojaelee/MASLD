#!/bin/bash
#SBATCH --job-name=qccontract
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay/logs/qccontract_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay/logs/qccontract_%j.err

set -eo pipefail
set +u
eval "$(micromamba shell hook -s bash)"
micromamba activate /gpfs/commons/home/jameslee/micromamba/envs/rnaseq
set -u

PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
SCRIPT_ROOT=${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay
: "${PLAN45_CANDIDATE_ID:?PLAN45_CANDIDATE_ID must name a new QC/unblinding contract}"
: "${PLAN45_STAGE_A_EXECUTION_CONTRACT_ROOT:?PLAN45_STAGE_A_EXECUTION_CONTRACT_ROOT is required}"
export MASLD_PROJECT_ROOT=${PROJECT_ROOT}
mkdir -p "${SCRIPT_ROOT}/logs"
cd "${SCRIPT_ROOT}"
python 33_build_stage_a_qc_unblinding_contract.py
python 34_validate_stage_a_qc_unblinding_contract.py
