#!/bin/bash
#SBATCH --job-name=relayreg
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay/logs/relayreg_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay/logs/relayreg_%j.err

set -eo pipefail
set +u
eval "$(micromamba shell hook -s bash)"
micromamba activate /gpfs/commons/home/jameslee/micromamba/envs/rnaseq
set -u

PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
SCRIPT_ROOT=${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay
: "${PLAN45_REGISTRY_CANDIDATE_ID:?PLAN45_REGISTRY_CANDIDATE_ID is required}"
: "${PLAN45_ORIENTATION_CANDIDATE_ID:?PLAN45_ORIENTATION_CANDIDATE_ID is required}"
export MASLD_PROJECT_ROOT=${PROJECT_ROOT}
mkdir -p "${SCRIPT_ROOT}/logs"

cd "${SCRIPT_ROOT}"
export PLAN45_CANDIDATE_ID=${PLAN45_REGISTRY_CANDIDATE_ID}
python 02_build_corrected_registry_v3.py
python 03_validate_registry_v3.py

export PLAN45_REGISTRY_ROOT=${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/${PLAN45_REGISTRY_CANDIDATE_ID}
export PLAN45_CANDIDATE_ID=${PLAN45_ORIENTATION_CANDIDATE_ID}
python 04_prepare_orientation_worklist.py
python 05_validate_orientation_worklist.py
