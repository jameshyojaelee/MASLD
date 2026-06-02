#!/bin/bash
#SBATCH --job-name=211_prog_driver_gen
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=results/stratified_causal/logs/211_%j.out
#SBATCH --error=results/stratified_causal/logs/211_%j.err

# Activate environment BEFORE set -u (activate-binutils references unbound ADDR2LINE)
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p results/stratified_causal/logs

echo "=== 211_progression_driver_genetics ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo ""

Rscript 211_progression_driver_genetics.R

echo ""
echo "End: $(date)"
