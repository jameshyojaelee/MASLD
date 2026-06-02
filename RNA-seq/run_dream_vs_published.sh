#!/bin/bash
#SBATCH --job-name=dream_vs_pub
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/dream_vs_pub_%j.out
#SBATCH --error=logs/dream_vs_pub_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript RNA-seq/206d_dream_vs_published.R
