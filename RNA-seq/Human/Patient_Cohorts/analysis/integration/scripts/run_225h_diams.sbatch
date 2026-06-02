#!/bin/bash
#SBATCH --job-name=225h_diams
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/225h_diams_%j.out
#SBATCH --error=logs/225h_diams_%j.err

set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"
LOG_DIR="$SCRIPTS_DIR/logs"
mkdir -p "$LOG_DIR"

echo "[$(date)] Starting 225h_diams_within_cohort.py"
echo "SLURM_JOB_ID=$SLURM_JOB_ID  CPUS=$SLURM_CPUS_PER_TASK"

micromamba run -n spatial python3 "$SCRIPTS_DIR/225h_diams_within_cohort.py"

echo "[$(date)] Done"
