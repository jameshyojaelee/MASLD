#!/bin/bash
#SBATCH --job-name=343c_fstage_alt
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343c_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343c_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

mkdir -p Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs

echo "[$(date)] starting 343c on $(hostname)"

micromamba run -n rapids_singlecell python -u \
    Analysis/SingleCell/scripts/343c_fstage_alt_inference_compare.py

echo "[$(date)] DONE"
