#!/bin/bash
#SBATCH --job-name=lineagefdr
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay/logs/lineagefdr_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay/logs/lineagefdr_%j.err

set -eo pipefail
set +u
eval "$(micromamba shell hook -s bash)"
micromamba activate /gpfs/commons/home/jameslee/micromamba/envs/rnaseq
set -u

PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
SCRIPT_ROOT=${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay
CANDIDATE_ID=source-independent-risk-state-relay-lineage-disease-reference-cv-fdr-correction-2026-08-10
mkdir -p "${SCRIPT_ROOT}/logs"
export MASLD_PROJECT_ROOT=${PROJECT_ROOT}
export PLAN45_DISEASE_CV_CORRECTION_ID=${CANDIDATE_ID}

python "${SCRIPT_ROOT}/20_correct_lineage_cv_multiplicity.py"
python "${SCRIPT_ROOT}/21_validate_lineage_cv_multiplicity.py"
