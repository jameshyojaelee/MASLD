#!/bin/bash
#SBATCH --job-name=209_sex_sub_fig
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/results/stratified_causal/logs/209_%j.out
#SBATCH --error=RNA-seq/results/stratified_causal/logs/209_%j.err

# Activate environment BEFORE set -u (activate-binutils references unbound ADDR2LINE)
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# Create log directory
mkdir -p RNA-seq/results/stratified_causal/logs

echo "=== Script 209: Sex-Subtype Synthesis Figures ==="
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

Rscript RNA-seq/209_sex_subtype_figures.R

echo "Finished: $(date)"
