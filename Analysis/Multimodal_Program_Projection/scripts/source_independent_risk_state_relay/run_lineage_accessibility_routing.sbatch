#!/usr/bin/env bash
#SBATCH --job-name=lineroute
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/lineroute_%j.out
#SBATCH --error=logs/lineroute_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPT_ROOT=${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay
RNASEQ_ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq

: "${PLAN45_CANDIDATE_ID:?PLAN45_CANDIDATE_ID is required}"
: "${PLAN45_ORIENTATION_ROOT:?PLAN45_ORIENTATION_ROOT is required}"
: "${PLAN45_LINEAGE_ROOT:?PLAN45_LINEAGE_ROOT is required}"
: "${PLAN45_ACCESSIBILITY_ROOT:?PLAN45_ACCESSIBILITY_ROOT is required}"

cd "${SCRIPT_ROOT}"
micromamba run -p "${RNASEQ_ENV}" Rscript 12_intersect_oriented_shared_signals.R
micromamba run -p "${RNASEQ_ENV}" python 13_validate_lineage_accessibility_routing.py
