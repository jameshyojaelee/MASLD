#!/bin/bash
#SBATCH --job-name=snapatac2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/peak_landscape_%j.out
#SBATCH --error=logs/peak_landscape_%j.err
set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate spatial
export PYTHONNOUSERSITE=1

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/ATAC/Human_Multiome
mkdir -p logs results/snapatac2

echo "==== Per-CT scATAC peak landscape annotation ===="
python scripts/31_export_peak_landscape.py
echo "==== Done ===="
