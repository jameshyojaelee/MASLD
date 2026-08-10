#!/bin/bash
#SBATCH --job-name=clccfix
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/logs/clccfix_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/logs/clccfix_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/logs
micromamba run -n rnaseq python Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge/13_prepare_clcc1_correction.py
export RISK_STATE_CANDIDATE_ID=chronic-state-risk-bridge-2026-08-09/clcc1_integrity_correction
micromamba run -n rnaseq Rscript Analysis/Multimodal_Program_Projection/scripts/risk_state_double_dissociation/05_run_clcc1.R
