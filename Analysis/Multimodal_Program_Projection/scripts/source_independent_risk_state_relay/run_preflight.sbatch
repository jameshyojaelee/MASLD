#!/bin/bash
#SBATCH --job-name=relayprep
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/source-independent-risk-state-relay-2026-08-10/logs/relayprep_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/source-independent-risk-state-relay-2026-08-10/logs/relayprep_%j.err

set -eo pipefail
set +u
eval "$(micromamba shell hook -s bash)"
micromamba activate /gpfs/commons/home/jameslee/micromamba/envs/rnaseq
set -u

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPT_ROOT=${ROOT}/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay
mkdir -p "${ROOT}/Analysis/Multimodal_Program_Projection/candidates/source-independent-risk-state-relay-2026-08-10/logs"
python "${SCRIPT_ROOT}/00_preflight_candidates.py"
python "${SCRIPT_ROOT}/01_validate_preflight.py"
