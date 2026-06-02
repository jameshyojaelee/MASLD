#!/bin/bash
#SBATCH --job-name=343o_heldout_cohort
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343o_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343o_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

mkdir -p Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs

echo "[$(date)] starting 343o on $(hostname)"

micromamba run -n rapids_singlecell python -u \
    Analysis/SingleCell/scripts/343o_heldout_cohort_replication.py

echo "[$(date)] DONE"
