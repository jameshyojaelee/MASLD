#!/bin/bash
#SBATCH --job-name=figS_transition_tf
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/figS_transition_tf_%j.out
#SBATCH --error=logs/figS_transition_tf_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p scripts/figures/logs

echo "=== figS_transition_tf_activity ==="
echo "Start: $(date)"
echo "Node: $(hostname)"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript scripts/figures/figS_transition_tf_activity.R

echo "End: $(date)"
