#!/bin/bash
#SBATCH --job-name=sex_v3_05
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=logs/M05_%j.out
#SBATCH --error=logs/M05_%j.err
# ---------------------------------------------------------------------------
# Module 05 — mashr fit on dream M2 betas/vcov.
# Wraps: scripts/sex_v3/05_mashr_fit.R
# Reads: intermediates/dream_M2_v3.rds  Writes: mashr_fit.rds + loadings/cov
# ---------------------------------------------------------------------------
set -eo pipefail

SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
# mashr is installed in the motifbreakr env, NOT rnaseq (fixed 2026-05-28)
micromamba activate motifbreakr

set -u
cd "${SCRIPT_DIR}"

echo "=== Module 05 start: $(date) ==="
echo "Host: $(hostname)  SLURM_JOB_ID: ${SLURM_JOB_ID:-NA}  CPUs: ${SLURM_CPUS_PER_TASK}"

Rscript 05_mashr_fit.R
if [ $? -ne 0 ]; then echo "FAILED: 05_mashr_fit.R"; exit 1; fi

echo "=== Module 05 end: $(date) ==="
