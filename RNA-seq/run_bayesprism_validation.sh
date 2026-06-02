#!/bin/bash
#SBATCH --job-name=liver_bayesprism_valid
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=06:00:00
#SBATCH --output=results/causal_inference/logs/bayesprism_validation_%j.out
#SBATCH --error=results/causal_inference/logs/bayesprism_validation_%j.err

# Script 22: BayesPrism Deconvolution Concordance & Attribution Validation
# Validates the 83% hepatocyte-intrinsic finding using BayesPrism (dual-method)
# Requires: merged_dge.rds, per-cohort BayesPrism/MuSiC proportions, dream_results.csv

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p results/causal_inference/logs
mkdir -p Human/Patient_Cohorts/analysis/integration/results/deconvolution/bayesprism
mkdir -p figures

echo "=== Script 22: BayesPrism Validation ==="
echo "Start: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-32}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript Human/Patient_Cohorts/analysis/integration/scripts/22_bayesprism_deconvolution.R

echo "=== Done: $(date) ==="
