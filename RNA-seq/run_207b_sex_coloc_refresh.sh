#!/bin/bash
#SBATCH --job-name=207b_sex_v3
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/results/stratified_causal/logs/207b_%j.out
#SBATCH --error=RNA-seq/results/stratified_causal/logs/207b_%j.err

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p RNA-seq/results/stratified_causal/logs

echo "=== Script 207b: Sex × COLOC refresh (v3 mashr sex_class) ==="
echo "Started: $(date)"
echo "Node: $(hostname)"

Rscript RNA-seq/207b_sex_coloc_refresh.R

echo "Finished: $(date)"
