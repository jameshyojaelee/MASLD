#!/bin/bash
#SBATCH --job-name=sex_v3_04
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=logs/M04_%j.out
#SBATCH --error=logs/M04_%j.err
# ---------------------------------------------------------------------------
# Module 04 — refit dream M2 (sex interaction) on -s 2 inputs.
# Wraps: scripts/sex_v3/04_dream_M2_extended.R
# Reads: intermediates/{sex_v3_input,sva_factors,age_mi}.rds  Writes: dream_M2_v3.rds
# ---------------------------------------------------------------------------
set -eo pipefail

SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u
cd "${SCRIPT_DIR}"

echo "=== Module 04 start: $(date) ==="
echo "Host: $(hostname)  SLURM_JOB_ID: ${SLURM_JOB_ID:-NA}  CPUs: ${SLURM_CPUS_PER_TASK}"

Rscript 04_dream_M2_extended.R
if [ $? -ne 0 ]; then echo "FAILED: 04_dream_M2_extended.R"; exit 1; fi

echo "=== Module 04 end: $(date) ==="
