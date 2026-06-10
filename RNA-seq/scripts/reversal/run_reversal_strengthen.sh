#!/bin/bash
#SBATCH --job-name=fgsea
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/reversal_strengthen_%j.out
#SBATCH --error=logs/reversal_strengthen_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/reversal
mkdir -p logs
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
echo "[$(date)] starting reversal gene-set strengthening on $(hostname)"
Rscript 02_reversal_geneset_strengthen.R
echo "[$(date)] done"
