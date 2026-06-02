#!/bin/bash
#SBATCH --job-name=343q_boot_stab
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --output=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343q_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343q_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

mkdir -p Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs
mkdir -p figures/supplementary/stage_ccc

echo "[$(date)] starting 343q on $(hostname)"
echo "SLURM_CPUS_PER_TASK = ${SLURM_CPUS_PER_TASK:-unset}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

micromamba run -n rapids_singlecell python -u \
    Analysis/SingleCell/scripts/343q_bootstrap_stability.py

echo "[$(date)] 343q finished"
