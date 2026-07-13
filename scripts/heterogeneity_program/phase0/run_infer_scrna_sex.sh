#!/bin/bash
#SBATCH --job-name=scanpy
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=200G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=scripts/heterogeneity_program/phase0/logs/scrna_sex_%j.out
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p scripts/heterogeneity_program/phase0/logs
source ~/.bashrc
micromamba activate spatial
python scripts/heterogeneity_program/phase0/infer_scrna_sex.py
