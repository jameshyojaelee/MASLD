#!/bin/bash
#SBATCH --job-name=218a_olink_triple
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/results/stratified_causal/logs/218a_%j.out
#SBATCH --error=RNA-seq/results/stratified_causal/logs/218a_%j.err

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p RNA-seq/results/stratified_causal/logs

echo "=== Script 218a: Sex × Olink × COLOC triple-concordance ==="
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-4}"

Rscript RNA-seq/218a_sex_olink_triple.R

echo "Finished: $(date)"
