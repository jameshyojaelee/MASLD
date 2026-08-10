#!/bin/bash
#SBATCH --job-name=identity
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/logs/identity_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/logs/identity_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/logs
micromamba run -n rnaseq python Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge/14_build_spatial_support.py
micromamba run -n rapids_singlecell python Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge/15_audit_program_identity.py
