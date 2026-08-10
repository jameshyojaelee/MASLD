#!/bin/bash
#SBATCH --job-name=stageatpl
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay/logs/stageatpl_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay/logs/stageatpl_%j.err

set -eo pipefail
set +u
eval "$(micromamba shell hook -s bash)"
micromamba activate /gpfs/commons/home/jameslee/micromamba/envs/rnaseq
set -u

PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
SCRIPT_ROOT=${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay
CANDIDATE_ID=source-independent-risk-state-relay-stage-a-template-v2-2026-08-10
mkdir -p "${SCRIPT_ROOT}/logs"
export MASLD_PROJECT_ROOT=${PROJECT_ROOT}
export PLAN45_STAGE_A_TEMPLATE_ID=${CANDIDATE_ID}

python "${SCRIPT_ROOT}/22_build_stage_a_template.py"
python "${SCRIPT_ROOT}/23_validate_stage_a_template.py"
