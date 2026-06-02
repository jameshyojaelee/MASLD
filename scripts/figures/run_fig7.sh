#!/bin/bash
#SBATCH --job-name=liver_fig7_multi_evidence
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/fig7_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/fig7_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "=== Generating Figure 7: Multi-Evidence Atlas ==="
echo "Started at $(date)"

Rscript scripts/figures/fig7_multi_evidence.R

echo "=== Done at $(date) ==="
