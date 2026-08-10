#!/bin/bash
#SBATCH --job-name=adjudicate
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/logs/adjudicate_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/logs/adjudicate_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
micromamba run -n rnaseq Rscript Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge/tests/test_state_axis.R
micromamba run -n rnaseq python Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge/16_adjudicate_plan43.py
micromamba run -n rnaseq python Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge/17_validate_plan43.py
