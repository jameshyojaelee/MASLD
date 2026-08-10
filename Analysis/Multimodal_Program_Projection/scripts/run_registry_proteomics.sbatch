#!/usr/bin/env bash
#SBATCH --job-name=proteomics
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/proteomics_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/proteomics_%j.err

set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
cd "$MASLD_PROJECT_ROOT"
"$R_BIN" Analysis/Multimodal_Program_Projection/scripts/01_freeze_program_registry.R
"$R_BIN" Analysis/Multimodal_Program_Projection/scripts/02_proteomics_projection.R
"$R_BIN" Analysis/Multimodal_Program_Projection/scripts/02b_histology_burden.R
