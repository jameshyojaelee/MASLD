#!/bin/bash
#SBATCH --job-name=hnf4a_zonation
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=logs/312e_zonation_%j.out
#SBATCH --error=logs/312e_zonation_%j.err

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate spatial

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell
mkdir -p logs results_gpu_v2/hepatocyte_subtypes/zonation_resolved

python scripts/312e_hnf4a_zonation_resolved.py
