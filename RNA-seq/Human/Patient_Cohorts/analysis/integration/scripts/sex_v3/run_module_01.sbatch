#!/bin/bash
#SBATCH --job-name=s_v3_M01
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=logs/M01_%j.out
#SBATCH --error=logs/M01_%j.err
# ---------------------------------------------------------------------------
# Module 01 — load + harmonize inputs for sex_v3.
# Wraps: scripts/sex_v3/01_load_data.R
# ---------------------------------------------------------------------------
set -e

SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "=== Module 01 start: $(date) ==="
echo "Host: $(hostname)  SLURM_JOB_ID: ${SLURM_JOB_ID:-NA}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}  Mem: 96G"
echo

Rscript 01_load_data.R

echo
echo "=== Module 01 end: $(date) ==="
