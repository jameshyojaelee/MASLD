#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --job-name=330_percell
#SBATCH --output=Analysis/SingleCell/logs/330_percell_%j.out
#SBATCH --error=Analysis/SingleCell/logs/330_percell_%j.err
set -eo pipefail; mkdir -p Analysis/SingleCell/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rapids_singlecell
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
python Analysis/SingleCell/scripts/330_percell_biological_signatures.py
