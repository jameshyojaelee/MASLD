#!/bin/bash
#SBATCH --job-name=218b_mendelian
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/results/stratified_causal/logs/218b_%j.out
#SBATCH --error=RNA-seq/results/stratified_causal/logs/218b_%j.err

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p RNA-seq/results/stratified_causal/logs

echo "=== Script 218b: Sex × Mendelian × COLOC overlap ==="
echo "Started: $(date)"
echo "Node: $(hostname)"

Rscript RNA-seq/218b_sex_mendelian_overlap.R

echo "Finished: $(date)"
