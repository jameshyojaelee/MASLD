#!/bin/bash
#SBATCH --job-name=lr_trajec
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/lr_traj_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/lr_traj_%j.err

set -euo pipefail
cd "${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
source ~/.bashrc 2>/dev/null || true

mkdir -p Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs

echo "[$(date)] Host: $(hostname)  CPUs: ${SLURM_CPUS_PER_TASK:-1}"

micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/352_lr_stage_trajectories.R

echo "[$(date)] Done."
