#!/bin/bash
#SBATCH --job-name=v6_8B
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/M07b2_%j.out
#SBATCH --error=logs/M07b2_%j.err
# ---------------------------------------------------------------------------
# Pillar 8B wrapper -- magnitude-aware LOCO concordance with bootstrap CIs.
# Reads intermediates/loco_v5/fold_{1..5}.rds; emits
# loco_concordance_per_class_v6.csv.
# ---------------------------------------------------------------------------
set -e

SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export LOCO_N_BOOT=${LOCO_N_BOOT:-1000}

echo "=== Pillar 8B start: $(date) ==="
echo "Host: $(hostname)  JID: ${SLURM_JOB_ID:-NA}  CPUs: ${SLURM_CPUS_PER_TASK}"
echo "N_BOOT=${LOCO_N_BOOT}"

Rscript 07b2_v6_loco_magnitude.R

echo "=== Pillar 8B end: $(date) ==="
