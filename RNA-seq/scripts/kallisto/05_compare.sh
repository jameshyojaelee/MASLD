#!/bin/bash
#SBATCH --job-name=B1_kallisto_compare
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/compare_%j.out
#SBATCH --error=logs/compare_%j.err

set -eo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

WT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
cd "$WT"
Rscript RNA-seq/scripts/kallisto/05_compare_kallisto_vs_star.R
