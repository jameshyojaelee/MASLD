#!/bin/bash
#SBATCH --job-name=quant_diag
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/quant_diag_%j.out
#SBATCH --error=logs/quant_diag_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript RNA-seq/206e_quantification_diagnostic.R
