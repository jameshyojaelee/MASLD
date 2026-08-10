#!/usr/bin/env bash
#SBATCH --job-name=riskstate
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/risk-state-double-dissociation-2026-08-09/logs/final_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/risk-state-double-dissociation-2026-08-09/logs/final_%j.err

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/risk_state_double_dissociation"
cd "$PROJECT_ROOT"
micromamba run -n rnaseq python "${SCRIPT_ROOT}/07_build_verdict.py"
micromamba run -n rnaseq python "${SCRIPT_ROOT}/08_validate_release.py"
