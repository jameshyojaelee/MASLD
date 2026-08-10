#!/bin/bash
#SBATCH --job-name=protein
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/logs/protein_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/logs/protein_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/logs
micromamba run -n rnaseq python Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge/09_build_gene_covariates.py
micromamba run -n rnaseq Rscript Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge/10_build_protein_adapters.R
micromamba run -n rnaseq Rscript Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge/11_run_protein_transportability.R
