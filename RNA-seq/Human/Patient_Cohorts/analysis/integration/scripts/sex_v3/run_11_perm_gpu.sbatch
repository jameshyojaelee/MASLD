#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --qos=nslab
#SBATCH --array=51-100%20
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --job-name=v6_P2_perm_gpu
#SBATCH --output=logs/P2_perm_gpu_%A_%a.out
#SBATCH --error=logs/P2_perm_gpu_%A_%a.err
set -e
SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export REP=${SLURM_ARRAY_TASK_ID}
Rscript 11_perm_fdr.R
