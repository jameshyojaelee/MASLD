#!/bin/bash
#SBATCH --job-name=net_295_indexes
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/net_295_indexes_%j.out
#SBATCH --error=RNA-seq/logs/net_295_indexes_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p RNA-seq/logs

export MAMBA_EXE='/gpfs/commons/home/jameslee/.local/bin/micromamba'
export MAMBA_ROOT_PREFIX='/gpfs/commons/home/jameslee/micromamba'
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate spatial

python3 RNA-seq/295_build_filter_indexes.py
