#!/bin/bash
#SBATCH --job-name=sex_v3_06
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=logs/M06_%j.out
#SBATCH --error=logs/M06_%j.err
# ---------------------------------------------------------------------------
# Module 06 — classify sex-dimorphic genes (v3) from mashr posteriors.
# Wraps: scripts/sex_v3/06_classify_v3.R
# Reads: intermediates/mashr_fit.rds + loadings/posteriors
# ---------------------------------------------------------------------------
set -eo pipefail

SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u
cd "${SCRIPT_DIR}"

echo "=== Module 06 start: $(date) ==="
echo "Host: $(hostname)  SLURM_JOB_ID: ${SLURM_JOB_ID:-NA}  CPUs: ${SLURM_CPUS_PER_TASK}"

Rscript 06_classify_v3.R
if [ $? -ne 0 ]; then echo "FAILED: 06_classify_v3.R"; exit 1; fi

echo "=== Module 06 end: $(date) ==="
