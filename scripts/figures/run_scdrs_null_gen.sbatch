#!/bin/bash
#SBATCH --job-name=scdrs_null_gen
#SBATCH --partition=cpu
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --output=scripts/figures/logs/scdrs_null_gen_%j.out
#SBATCH --error=scripts/figures/logs/scdrs_null_gen_%j.err

set -eo pipefail
mkdir -p scripts/figures/logs

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate spatial

python scripts/figures/figS_scdrs_null_gen.py
