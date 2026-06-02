#!/bin/bash
#SBATCH --job-name=pub_vs_perstudy
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=logs/pub_vs_perstudy_%j.out
#SBATCH --error=logs/pub_vs_perstudy_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript RNA-seq/206_published_vs_perstudy.R
