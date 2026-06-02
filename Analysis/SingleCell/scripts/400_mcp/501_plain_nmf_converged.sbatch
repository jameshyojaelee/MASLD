#!/bin/bash
#SBATCH --job-name=alpha1_plain_nmf
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=400G
#SBATCH --cpus-per-task=16
#SBATCH --time=72:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/501_plain_nmf_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/501_plain_nmf_%j.err
set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp
source "$HOME/.bashrc"
micromamba activate rapids_singlecell

# Thread knobs so BLAS can use the full allocation within each NMF fit
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

# Ensure fresh runs can be tuned from the sbatch env
: "${PLAIN_N_REP:=50}"
: "${PLAIN_K_LIST:=10,16}"
: "${PLAIN_MAX_ITER:=1000}"
: "${PLAIN_TOL:=1e-4}"
: "${PLAIN_N_CELLS:=0}"
: "${PLAIN_CONSENSUS_N_CELLS:=30000}"
export PLAIN_N_REP PLAIN_K_LIST PLAIN_MAX_ITER PLAIN_TOL PLAIN_N_CELLS PLAIN_CONSENSUS_N_CELLS

python -u 501_plain_nmf_converged.py
