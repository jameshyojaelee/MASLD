#!/bin/bash
#SBATCH --job-name=218c_sex_gsea
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/results/stratified_causal/logs/218c_%j.out
#SBATCH --error=RNA-seq/results/stratified_causal/logs/218c_%j.err

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p RNA-seq/results/stratified_causal/logs

echo "=== Script 218c: Sex × COLOC pathway fgsea ==="
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-8}"

Rscript RNA-seq/218c_sex_coloc_pathway_gsea.R

echo "Finished: $(date)"
