#!/bin/bash
#SBATCH --job-name=extract_sdc2_C4
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --mem=160G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/extract_sdc2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/extract_sdc2_%j.err

set -euo pipefail
mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

source /gpfs/commons/home/jameslee/.bashrc 2>/dev/null || true
eval "$(micromamba shell hook --shell bash)"
micromamba activate spatial

python Analysis/SingleCell/scripts/400_mcp/extract_sdc2_per_cell.py
