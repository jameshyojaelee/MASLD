#!/bin/bash
#SBATCH --job-name=212_progression_figs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/212_progression_figs_%j.out
#SBATCH --error=logs/212_progression_figs_%j.err

# Activate environment BEFORE set -u (activate-binutils references unbound ADDR2LINE)
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -euo pipefail

echo "=== Script 212: Progression Figures + Atlas Integration ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo ""

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

cd "$MASLD_PROJECT_ROOT/RNA-seq"
mkdir -p logs

Rscript 212_progression_figures.R

echo ""
echo "Completed: $(date)"
