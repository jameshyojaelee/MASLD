#!/bin/bash
#SBATCH --job-name=155_prog_fig
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=../results/prognosis/155_figures_%j.out
#SBATCH --error=../results/prognosis/155_figures_%j.err

set -euo pipefail

echo "=== 155_prognosis_figures ==="
echo "Host: $(hostname)"
echo "Date: $(date)"
echo "SLURM Job ID: ${SLURM_JOB_ID}"

mkdir -p ../results/prognosis

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

Rscript 155_prognosis_figures.R

echo "Done: $(date)"
