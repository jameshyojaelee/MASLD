#!/bin/bash
#SBATCH --job-name=343d_fstage_val
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343d_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343d_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

mkdir -p Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs
mkdir -p figures/supplementary/stage_ccc

echo "[$(date)] starting 343d on $(hostname)"

micromamba run -n rapids_singlecell python -u \
    Analysis/SingleCell/scripts/343d_fstage_external_validation.py

echo "[$(date)] DONE"
