#!/bin/bash
#SBATCH --job-name=s_v3_M02
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=logs/M02_%j.out
#SBATCH --error=logs/M02_%j.err
# ---------------------------------------------------------------------------
# Module 02 — SVA hidden factor estimation for sex_v3.
# Wraps: scripts/sex_v3/02_sva_factors.R
# ---------------------------------------------------------------------------
set -e

SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "=== Module 02 start: $(date) ==="
echo "Host: $(hostname)  SLURM_JOB_ID: ${SLURM_JOB_ID:-NA}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}  Mem: 96G"
echo

Rscript 02_sva_factors.R

echo
echo "=== Module 02 end: $(date) ==="
