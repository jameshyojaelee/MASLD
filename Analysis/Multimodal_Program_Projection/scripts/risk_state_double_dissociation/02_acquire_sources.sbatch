#!/usr/bin/env bash
#SBATCH --job-name=source
#SBATCH --partition=io
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/risk-state-double-dissociation-2026-08-09/logs/source_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/risk-state-double-dissociation-2026-08-09/logs/source_%j.err

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "$PROJECT_ROOT"
python3 Analysis/Multimodal_Program_Projection/scripts/risk_state_double_dissociation/02_acquire_sources.py

