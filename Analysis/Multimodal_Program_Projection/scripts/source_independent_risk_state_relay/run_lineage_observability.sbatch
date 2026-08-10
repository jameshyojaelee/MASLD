#!/bin/bash
#SBATCH --job-name=lineageobs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay/logs/lineageobs_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay/logs/lineageobs_%j.err

set -eo pipefail
set +u
eval "$(micromamba shell hook -s bash)"
micromamba activate /gpfs/commons/home/jameslee/micromamba/envs/rnaseq
set -u

PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
SCRIPT_ROOT=${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay
CANDIDATE_ID=source-independent-risk-state-relay-lineage-observability-2026-08-10
mkdir -p "${SCRIPT_ROOT}/logs"
export MASLD_PROJECT_ROOT=${PROJECT_ROOT}
export PLAN45_CANDIDATE_ID=${CANDIDATE_ID}
export PLAN45_LINEAGE_CANDIDATE_ID=${CANDIDATE_ID}

Rscript "${SCRIPT_ROOT}/08_build_donor_lineage_observability.R"
python "${SCRIPT_ROOT}/09_validate_lineage_observability.py"
