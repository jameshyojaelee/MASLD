#!/bin/bash
#SBATCH --job-name=script63_ordinal
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/script63_%j.out
#SBATCH --error=logs/script63_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate rnaseq

echo "=== Script 63: Ordinal Models (per-fold genes) ==="
echo "Started: $(date)"
Rscript 63_ordinal_models.R 2>&1
echo "Completed: $(date)"
