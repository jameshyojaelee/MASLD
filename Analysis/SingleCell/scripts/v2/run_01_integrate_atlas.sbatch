#!/bin/bash
#SBATCH --job-name=v2_integrate_atlas
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=24:00:00
#SBATCH --output=Analysis/SingleCell/scripts/v2/logs/01_integrate_%j.log
#SBATCH --error=Analysis/SingleCell/scripts/v2/logs/01_integrate_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

mkdir -p Analysis/SingleCell/scripts/v2/logs

echo "[$(date)] Job ${SLURM_JOB_ID} starting on $(hostname)"
echo "[$(date)] CPUs: ${SLURM_CPUS_PER_TASK}"

export PYTHONHASHSEED=42
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

micromamba run -n rapids_singlecell python \
  Analysis/SingleCell/scripts/v2/01_integrate_atlas_v2.py

echo "[$(date)] DONE"
