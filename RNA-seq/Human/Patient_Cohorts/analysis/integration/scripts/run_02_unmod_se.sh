#!/bin/bash
#SBATCH --job-name=limma-voom
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/limma_voom_%j.out
#SBATCH --error=logs/limma_voom_%j.err

set -eo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts

echo "=== Script 02: Per-study DE with unmoderated SE ==="
echo "Started: $(date)"
echo "Node: $(hostname)"

Rscript 02_per_study_de.R

echo "=== Script 02 complete: $(date) ==="
