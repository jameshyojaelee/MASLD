#!/bin/bash
#SBATCH --job-name=plasma_etiology_225e
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/225e_plasma_etiology_%j.out
#SBATCH --error=logs/225e_plasma_etiology_%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
mkdir -p logs

echo "=== 225e plasma etiology sweep ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"

micromamba run -n spatial python 225e_plasma_sweep_etiology.py

echo "End: $(date)"
