#!/bin/bash
#SBATCH --job-name=schloproj
#SBATCH --partition=cpu
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/logs/schlo_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/logs/schlo_%j.err

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/public_functional"
cd "${PROJECT_ROOT}"
micromamba run -n rapids_singlecell python "${SCRIPT_ROOT}/05_build_schlo_pseudobulk.py"
micromamba run -n rnaseq Rscript "${SCRIPT_ROOT}/06_run_schlo_programs.R"
