#!/bin/bash
#SBATCH --job-name=343s_dub_thr
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --output=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343s_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343s_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

mkdir -p Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs

echo "[$(date)] starting 343s on $(hostname)"

micromamba run -n rapids_singlecell python -u \
    Analysis/SingleCell/scripts/343s_dubious_threshold_sensitivity.py

echo "[$(date)] DONE"
