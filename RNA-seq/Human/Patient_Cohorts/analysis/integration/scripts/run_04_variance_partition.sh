#!/usr/bin/env bash
#SBATCH --job-name=varpart_04
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/04_variance_partition_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/04_variance_partition_%j.err

# R1-C13: Re-run variance partition with group_binary instead of condition
set -euo pipefail

SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "$LOG_DIR"

echo "[$(date)] Starting variance partition (R1-C13 fix: group_binary)"
echo "SLURM Job ID: ${SLURM_JOB_ID:-local}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-1}"

micromamba run -n rnaseq Rscript "${SCRIPT_DIR}/04_variance_partition.R"

echo "[$(date)] Done"
