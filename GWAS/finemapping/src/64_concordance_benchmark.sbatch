#!/bin/bash
#SBATCH --job-name=borzoi
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/64_concordance_%j.log
set -euo pipefail
export PYTHONHASHSEED=0 PYTHONNOUSERSITE=1
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/src
python 64_concordance_benchmark.py
