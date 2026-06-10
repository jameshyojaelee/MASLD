#!/bin/bash
#SBATCH --job-name=limma
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/reversal_de_%j.out
#SBATCH --error=logs/reversal_de_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/reversal
mkdir -p logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

echo "[$(date)] starting paired reversal DE on $(hostname)"
Rscript 01_paired_reversal_de.R
echo "[$(date)] done"
