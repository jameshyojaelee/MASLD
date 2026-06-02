#!/bin/bash
#SBATCH --job-name=pb_bulk_concordance
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=logs/36b_concordance_%j.log

set -euo pipefail

echo "=== Pseudobulk-Bulk Concordance (Phase 1 baseline) ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"

eval "$(micromamba shell hook --shell=bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p logs

Rscript 36b_pseudobulk_bulk_concordance.R

echo "End: $(date)"
