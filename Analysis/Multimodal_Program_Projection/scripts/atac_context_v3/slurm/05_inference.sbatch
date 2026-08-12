#!/bin/bash
#SBATCH --job-name=atacinf
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_inference_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_inference_%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PACKAGE=${ROOT}/Analysis/Multimodal_Program_Projection/scripts/atac_context_v3
MAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
eval "$(${MAMBA} shell hook --shell bash)"
micromamba activate rnaseq
Rscript "${PACKAGE}/04_run_da.R"
Rscript "${PACKAGE}/05_run_programs.R"
Rscript "${PACKAGE}/10b_condition_blindness_audit.R"
Rscript "${PACKAGE}/12_render_non_genetic.R"
