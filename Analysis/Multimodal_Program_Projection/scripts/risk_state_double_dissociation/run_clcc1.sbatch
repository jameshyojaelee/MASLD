#!/usr/bin/env bash
#SBATCH --job-name=clcc1
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/risk-state-double-dissociation-2026-08-09/logs/clcc1_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/risk-state-double-dissociation-2026-08-09/logs/clcc1_%j.err

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/risk_state_double_dissociation"
cd "$PROJECT_ROOT"
python3 "${SCRIPT_ROOT}/04_prepare_huh7_covariates.py"
micromamba run -n rnaseq Rscript "${SCRIPT_ROOT}/05_run_clcc1.R"
