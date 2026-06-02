#!/bin/bash
#SBATCH --job-name=switch_class
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/results/stratified_causal/logs/switch_classification_%j.log
#SBATCH --error=RNA-seq/results/stratified_causal/logs/switch_classification_%j.log

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p RNA-seq/results/stratified_causal/logs

echo "=== Switch gene classification (Gap #14) ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start:  $(date)"
echo "Node:   $(hostname)"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

Rscript RNA-seq/218_switch_gene_classification.R

echo "End: $(date)"
