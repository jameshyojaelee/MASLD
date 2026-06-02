#!/bin/bash
#SBATCH --job-name=ccc_fig_v2
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/353_heatmap_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/353_heatmap_%j.err

set -euo pipefail
cd "${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
source ~/.bashrc 2>/dev/null || true

mkdir -p Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs

echo "[$(date)] starting 353_enhanced_trajectory_heatmap.R on $(hostname)"
micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/353_enhanced_trajectory_heatmap.R
echo "[$(date)] done"
