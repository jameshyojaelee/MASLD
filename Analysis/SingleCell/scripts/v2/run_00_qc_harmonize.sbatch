#!/bin/bash
#SBATCH --job-name=v2_qc_harmonize
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=12:00:00
#SBATCH --output=Analysis/SingleCell/scripts/v2/logs/00_qc_harmonize_%j.log
#SBATCH --error=Analysis/SingleCell/scripts/v2/logs/00_qc_harmonize_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

mkdir -p Analysis/SingleCell/scripts/v2/logs

echo "[$(date)] Job ${SLURM_JOB_ID} starting on $(hostname)"
echo "[$(date)] CPUs: ${SLURM_CPUS_PER_TASK}"

export PYTHONHASHSEED=42
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

micromamba run -n rapids_singlecell python \
  Analysis/SingleCell/scripts/v2/00_qc_harmonize_per_dataset.py

echo "[$(date)] DONE"
