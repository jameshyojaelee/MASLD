#!/bin/bash
#SBATCH --job-name=regen_panels
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/regen_panels_%j.out
#SBATCH --error=logs/regen_panels_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p scripts/figures/logs

eval "$(micromamba shell hook --shell=bash)"
micromamba activate rapids_singlecell
export MPLBACKEND=Agg

python scripts/figures/regen_panels_1_5_10.py
