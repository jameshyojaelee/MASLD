#!/bin/bash
#SBATCH --job-name=343r_anchor_jk
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --output=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343r_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343r_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

mkdir -p Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs
mkdir -p figures/supplementary/stage_ccc

echo "[$(date)] starting 343r on $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"

micromamba run -n rapids_singlecell python -u \
    Analysis/SingleCell/scripts/343r_anchor_jackknife.py

echo "[$(date)] DONE"
