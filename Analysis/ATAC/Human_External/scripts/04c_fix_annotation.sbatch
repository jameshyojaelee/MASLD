#!/bin/bash
#SBATCH --job-name=snapatac2
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/ATAC/Human_External/logs/fixannot_%j.out
#SBATCH --error=Analysis/ATAC/Human_External/logs/fixannot_%j.err
set -eo pipefail
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
echo "[$(date)] 04c fix annotation on $(hostname)"
python Analysis/ATAC/Human_External/scripts/04c_fix_annotation.py
echo "[$(date)] done"
