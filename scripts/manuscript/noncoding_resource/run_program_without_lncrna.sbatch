#!/bin/bash
#SBATCH --job-name=nclnprog
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/nclnprog_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/nclnprog_%j.err

set -euo pipefail
umask 027

: "${NONCODING_CELL_TYPE:?Set NONCODING_CELL_TYPE}"
: "${NONCODING_PROGRAM_OUTPUT:?Set NONCODING_PROGRAM_OUTPUT}"

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export PATH="/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin:${PATH}"

/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python \
  scripts/manuscript/noncoding_resource/rescore_programs_without_lncrna.py \
  --cell-type "${NONCODING_CELL_TYPE}" \
  --output "${NONCODING_PROGRAM_OUTPUT}"
