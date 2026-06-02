#!/bin/bash
#SBATCH --job-name=v6_P3
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=logs/P3_%A_%a.out
#SBATCH --error=logs/P3_%A_%a.err
#SBATCH --array=1-60%12
# Pillar 3 — selection-bias simulation. 120 cells split: 1-60 on cpu, 61-120
# on gpu (see run_12_selection_bias_gpu.sbatch). TID is global so per-cell RDS
# paths collide cleanly into the same intermediates/selection_bias_v6 dir.
set -e
SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
echo "=== P3 task ${SLURM_ARRAY_TASK_ID} (cpu) ==="
echo "Host: $(hostname)  date: $(date)"
Rscript 12_selection_bias_sim.R
echo "=== done: $(date) ==="
