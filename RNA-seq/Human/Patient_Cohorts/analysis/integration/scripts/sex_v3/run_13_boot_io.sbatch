#!/bin/bash
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --array=41-50%10
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --job-name=v6_P4_boot_io
#SBATCH --output=logs/P4_boot_io_%A_%a.out
#SBATCH --error=logs/P4_boot_io_%A_%a.err
set -e
SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export REP=${SLURM_ARRAY_TASK_ID}
Rscript 13_bootstrap_dream.R
