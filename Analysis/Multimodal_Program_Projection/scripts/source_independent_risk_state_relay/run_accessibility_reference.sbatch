#!/usr/bin/env bash
#SBATCH --job-name=accessref
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/accessref_%j.out
#SBATCH --error=logs/accessref_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPT_ROOT=${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay

cd "${SCRIPT_ROOT}"
export PLAN45_CANDIDATE_ID=${PLAN45_CANDIDATE_ID:-source-independent-risk-state-relay-accessibility-reference-2026-08-10}

python3 10_build_accessibility_reference.py
python3 11_validate_accessibility_reference.py
