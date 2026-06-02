#!/bin/bash
#SBATCH --job-name=hoang_invest
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/hoang_invest_%j.out
#SBATCH --error=logs/hoang_invest_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript RNA-seq/206b_hoang_method_investigation.R
