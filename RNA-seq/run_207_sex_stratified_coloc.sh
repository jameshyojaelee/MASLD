#!/bin/bash
#SBATCH --job-name=207_sex_coloc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/results/stratified_causal/logs/207_%j.out
#SBATCH --error=RNA-seq/results/stratified_causal/logs/207_%j.err

# Activate environment BEFORE set -u (activate-binutils references unbound ADDR2LINE)
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# Create log directory
mkdir -p RNA-seq/results/stratified_causal/logs

echo "=== Script 207: Sex-Stratified COLOC Enrichment ==="
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

Rscript RNA-seq/207_sex_stratified_coloc.R

echo "Finished: $(date)"
