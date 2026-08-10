#!/usr/bin/env bash
#SBATCH --job-name=guidegate
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=logs/guidegate_%j.out
#SBATCH --error=logs/guidegate_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPT_ROOT=${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay

: "${PLAN45_CANDIDATE_ID:?PLAN45_CANDIDATE_ID is required}"
: "${PLAN45_ROUTING_ROOT:?PLAN45_ROUTING_ROOT is required}"

cd "${SCRIPT_ROOT}"
python3 14_prepare_guide_design_worklist.py
python3 15_validate_guide_design_worklist.py
