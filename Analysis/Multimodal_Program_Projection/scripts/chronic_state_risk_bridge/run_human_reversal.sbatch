#!/bin/bash
#SBATCH --job-name=humanrev
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/logs/humanrev_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/logs/humanrev_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/logs
micromamba run -n rnaseq python Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge/06_build_human_adapters.py
micromamba run -n rnaseq Rscript Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge/07_run_human_reversal.R
mkdir -p Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/execution
export PLAN43_EXACT_NAS_BINARY="$PWD/Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/execution/exact_nas_permutation"
g++ -O3 -std=c++17 Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge/exact_nas_permutation.cpp -o "$PLAN43_EXACT_NAS_BINARY"
micromamba run -n rnaseq Rscript Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge/08_finalize_human_reversal.R
