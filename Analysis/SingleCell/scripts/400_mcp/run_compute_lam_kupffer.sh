#!/bin/bash
#SBATCH --job-name=lam_kc_score
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/lam_kc_score_%j.out
#SBATCH --error=logs/lam_kc_score_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

mkdir -p Analysis/SingleCell/scripts/400_mcp/logs

source /gpfs/commons/home/jameslee/.bashrc
micromamba activate spatial

python Analysis/SingleCell/scripts/400_mcp/compute_lam_kupffer_scores.py
