#!/bin/bash
#SBATCH --job-name=panel_benchmark
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/panel_benchmark_%j.log
#SBATCH --error=logs/panel_benchmark_%j.log

set -euo pipefail

echo "=== Published Panel Benchmark ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $(hostname)"
echo "Start:  $(date)"

# Ensure log directory exists
mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs

# Activate environment
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$MASLD_PROJECT_ROOT"

Rscript scripts/figures/figS_published_panel_benchmark.R

echo "Done: $(date)"
