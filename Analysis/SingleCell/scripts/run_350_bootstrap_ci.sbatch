#!/bin/bash
#SBATCH --job-name=ccc_bootci
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/350_bootstrap_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/350_bootstrap_%j.err

set -euo pipefail
cd "${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
source ~/.bashrc 2>/dev/null || true

mkdir -p Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs

echo "[$(date)] starting 350_stage_ccc_bootstrap.R on $(hostname)"
echo "[$(date)] SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-?}"

micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/350_stage_ccc_bootstrap.R

echo "[$(date)] finished"
