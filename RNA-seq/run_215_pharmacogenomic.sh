#!/bin/bash
#SBATCH --job-name=s215_pharma
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/stratified_causal/logs/s215_pharmacogenomic_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/stratified_causal/logs/s215_pharmacogenomic_%j.err

# Activate environment BEFORE set -u (activate-binutils references unbound ADDR2LINE)
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -euo pipefail

echo "=== Script 215: Pharmacogenomic Target Mapping ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo ""

# Create log directory
mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/stratified_causal/logs

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

cd "${MASLD_PROJECT_ROOT}/RNA-seq"
Rscript 215_pharmacogenomic_targets.R

echo ""
echo "Completed: $(date)"
