#!/bin/bash
#SBATCH --job-name=s217_strat_atlas
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/217_stratified_causal_atlas_%j.out
#SBATCH --error=logs/217_stratified_causal_atlas_%j.err

# Activate environment BEFORE set -u (activate-binutils references unbound ADDR2LINE)
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -euo pipefail

echo "=== Script 217: Stratified Causal Atlas Integration ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo ""

# Set project root
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

cd "$MASLD_PROJECT_ROOT/RNA-seq"

# Ensure logs directory exists
mkdir -p logs

# Ensure output directory exists
mkdir -p results/stratified_causal

Rscript 217_stratified_causal_atlas.R

echo ""
echo "End: $(date)"
echo "=== Done ==="
