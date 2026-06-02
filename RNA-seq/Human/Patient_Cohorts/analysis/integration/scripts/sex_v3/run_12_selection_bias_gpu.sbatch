#!/bin/bash
#SBATCH --job-name=v6_P3g
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=logs/P3_gpu_%A_%a.out
#SBATCH --error=logs/P3_gpu_%A_%a.err
#SBATCH --array=61-120%12
set -e
SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
echo "=== P3 task ${SLURM_ARRAY_TASK_ID} (gpu) ==="
echo "Host: $(hostname)  date: $(date)"
Rscript 12_selection_bias_sim.R
echo "=== done: $(date) ==="
