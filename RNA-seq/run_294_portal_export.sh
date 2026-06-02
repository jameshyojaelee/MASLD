#!/bin/bash
#SBATCH --job-name=net_294_portal
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/net_294_portal_%j.out
#SBATCH --error=logs/net_294_portal_%j.err

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

echo "=== 294_export_portal_v2.py ==="
echo "Job ID: ${SLURM_JOB_ID:-interactive}"
echo "Start: $(date)"

python 294_export_portal_v2.py

echo "End: $(date)"
echo "=== done ==="
