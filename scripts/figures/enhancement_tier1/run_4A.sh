#!/bin/bash
#SBATCH --job-name=figures
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/fig4A_%j.out
#SBATCH --error=logs/fig4A_%j.err
set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/enhancement_tier1
mkdir -p logs
eval "$(micromamba shell hook --shell bash)"; micromamba activate rnaseq; set -u
echo "[$(date)] 4A orthogonality panels"; Rscript make_4A_orthogonality.R; echo "[$(date)] done"
