#!/bin/bash
#SBATCH --job-name=343m_aug_fstage
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343m_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343m_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

mkdir -p Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs

echo "[$(date)] starting 343m on $(hostname)"

micromamba run -n rapids_singlecell python -u \
    Analysis/SingleCell/scripts/343m_augmented_fstage.py

echo "[$(date)] DONE"
