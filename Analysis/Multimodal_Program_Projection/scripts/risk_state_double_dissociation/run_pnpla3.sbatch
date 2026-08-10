#!/usr/bin/env bash
#SBATCH --job-name=pnpla3
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/risk-state-double-dissociation-2026-08-09/logs/pnpla3_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/risk-state-double-dissociation-2026-08-09/logs/pnpla3_%j.err

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/risk_state_double_dissociation"
cd "$PROJECT_ROOT"
micromamba run -n rnaseq Rscript "${SCRIPT_ROOT}/06_run_pnpla3.R"

