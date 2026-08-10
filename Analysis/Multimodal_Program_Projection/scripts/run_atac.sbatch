#!/usr/bin/env bash
#SBATCH --job-name=ATAC
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/ATAC_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/ATAC_%j.err

set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
cd "$MASLD_PROJECT_ROOT"
"$R_BIN" Analysis/Multimodal_Program_Projection/scripts/03_atac_projection.R

