#!/bin/bash
#SBATCH --job-name=B4_uq_sens
#SBATCH --partition=bigmem
#SBATCH --mem=500G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=logs/B4/uq_sens_%j.out
#SBATCH --error=logs/B4/uq_sens_%j.err

set -eo pipefail
source /gpfs/commons/home/jameslee/.bashrc
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export R_PARALLEL_SEED=42

echo "Host: $(hostname); CPUS=${SLURM_CPUS_PER_TASK}"
Rscript RNA-seq/scripts/14_5_uq_normalization_sensitivity.R
echo "Done."
