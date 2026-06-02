#!/bin/bash
#SBATCH --job-name=v6_P5
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=logs/P5_%A_%a.out
#SBATCH --error=logs/P5_%A_%a.err
#SBATCH --array=1-100%16
# Pillar 5 — 50-seed calibration. TID in [1..200] = 50 seeds × 4 effect_mag.
# Split: tasks 1-100 on cpu, 101-200 on gpu (see _gpu.sbatch).
# Decode TID -> (SEED, EFFECT_MAG):
#   SEED       = ((TID - 1) / 4) + 1    (1..50)
#   EFF_IDX    = (TID - 1) %% 4          (0..3)
set -e
SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

TID=${SLURM_ARRAY_TASK_ID}
SEED=$(( ((TID - 1) / 4) + 1 ))
EFF_IDX=$(( (TID - 1) % 4 ))
case $EFF_IDX in
  0) EFFECT_MAG=0.30 ;;
  1) EFFECT_MAG=0.50 ;;
  2) EFFECT_MAG=0.80 ;;
  3) EFFECT_MAG=1.20 ;;
esac
export SEED EFFECT_MAG
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
echo "=== P5 task ${TID}: SEED=${SEED} EFFECT_MAG=${EFFECT_MAG} (cpu) ==="
echo "Host: $(hostname)  date: $(date)"
Rscript 14_calibration_50seed.R
echo "=== done: $(date) ==="
