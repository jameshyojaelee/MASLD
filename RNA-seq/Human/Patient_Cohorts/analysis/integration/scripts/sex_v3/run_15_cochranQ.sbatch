#!/bin/bash
#SBATCH --job-name=v6_P6
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/P6_%A_%a.out
#SBATCH --error=logs/P6_%A_%a.err
#SBATCH --array=1-5%5
set -e
SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
echo "=== P6 task ${SLURM_ARRAY_TASK_ID} (cpu, per-cohort fit) ==="
echo "Host: $(hostname)  date: $(date)"
Rscript 15_cochranQ_per_cohort.R
echo "=== done: $(date) ==="
