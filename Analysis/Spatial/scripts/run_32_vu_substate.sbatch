#!/bin/bash
#SBATCH --job-name=vu_substate
#SBATCH --output=logs/vu_substate_%j.log
#SBATCH --error=logs/vu_substate_%j.log
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -e
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate spatial

echo "=== Vu et al. spot-level F3a/F3b ==="
echo "Started: $(date)"

python scripts/32_vu_substate_f3a_f3b.py

echo "Finished: $(date)"
