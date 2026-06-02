#!/bin/bash
#SBATCH --job-name=347b_dissoc
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/347b_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/347b_%j.err

set -euo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$BASE"
mkdir -p Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs

export MASLD_PROJECT_ROOT="$BASE"

echo "[$(date)] Starting 347b dissociation sensitivity test"
micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/347b_gse136103_holdout_sensitivity.R
echo "[$(date)] Done"
