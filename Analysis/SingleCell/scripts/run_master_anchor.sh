#!/bin/bash
#SBATCH --job-name=anchor
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/master_anchor_%j.out
#SBATCH --error=logs/master_anchor_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts
mkdir -p logs
eval "$(micromamba shell hook --shell bash)"
micromamba activate spatial
set -u
echo "[$(date)] Tier 1C master anchor on $(hostname)"
python WS4_01_master_anchor.py
echo "[$(date)] done"
