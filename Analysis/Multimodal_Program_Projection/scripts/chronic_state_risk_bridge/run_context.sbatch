#!/bin/bash
#SBATCH --job-name=context
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/logs/context_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/logs/context_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/logs
micromamba run -n rnaseq Rscript Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge/12_run_context_topology.R
