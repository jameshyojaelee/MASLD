#!/bin/bash
#SBATCH --job-name=susieR
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --array=1-40
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/pip_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/pip_%A_%a.err
set -eo pipefail
set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate finemapping
set -u
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
Rscript /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/tmp_diag/pip_paired.R
