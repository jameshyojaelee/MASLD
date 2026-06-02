#!/bin/bash
#SBATCH --job-name=sex-permutation
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --mem=500G
#SBATCH --cpus-per-task=16
#SBATCH --time=72:00:00
#SBATCH --output=logs/sex_perm_%j.out
#SBATCH --error=logs/sex_perm_%j.err

set -eo pipefail
source /gpfs/commons/home/jameslee/.bashrc
micromamba activate rnaseq
set -u

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "${PROJ}"

export MASLD_PROJECT_ROOT="${PROJ}"
export R_PARALLEL_SEED=42
export N_PERM=${N_PERM:-10000}

echo "Host: $(hostname); CPUS=${SLURM_CPUS_PER_TASK}; N_PERM=${N_PERM}; START=$(date)"
Rscript "${PROJ}/RNA-seq/scripts/14_4_sex_permutation_null.R"
echo "Done. END=$(date)"
