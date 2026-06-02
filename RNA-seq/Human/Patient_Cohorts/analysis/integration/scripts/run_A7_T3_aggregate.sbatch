#!/bin/bash
#SBATCH --job-name=A7_T3_agg
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/A7_T3_agg_%j.out
#SBATCH --error=logs/A7_T3_agg_%j.err

set -eo pipefail
mkdir -p logs

echo "=== A7 T3 aggregate ==="
echo "Started: $(date)"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript 14.5a_T3_aggregate.R

echo "Done: $(date)"
