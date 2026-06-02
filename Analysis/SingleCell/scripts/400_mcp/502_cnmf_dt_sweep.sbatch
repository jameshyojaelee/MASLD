#!/bin/bash
#SBATCH --job-name=alpha2_cnmf_dt
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=750G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/502_cnmf_dt_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/502_cnmf_dt_%j.err
set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp
source "$HOME/.bashrc"
micromamba activate rapids_singlecell

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

: "${CNMF_DT_LIST:=0.022,0.03,0.05}"
export CNMF_DT_LIST

python -u 502_cnmf_dt_sweep.py
