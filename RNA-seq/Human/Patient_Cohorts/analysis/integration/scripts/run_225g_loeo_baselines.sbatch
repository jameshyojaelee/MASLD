#!/bin/bash
#SBATCH --job-name=plasma_loeo_225g
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/225g_plasma_loeo_%j.out
#SBATCH --error=logs/225g_plasma_loeo_%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
mkdir -p logs

echo "=== 225g plasma LOEO, random baselines, NFASC+GDF15 benchmark ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"

micromamba run -n spatial python 225g_plasma_loeo_baselines.py

echo "End: $(date)"
