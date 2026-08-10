#!/usr/bin/env bash
#SBATCH --job-name=Rscript
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/validate_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/validate_%j.err

set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
cd "$MASLD_PROJECT_ROOT"
"$R_BIN" Analysis/Multimodal_Program_Projection/scripts/06_validate_outputs.R

