#!/bin/bash
#SBATCH --job-name=liver_drug_repurpose
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=results/causal_inference/logs/drug_repurposing_%j.out
#SBATCH --error=results/causal_inference/logs/drug_repurposing_%j.err

# Strategy 11: Drug Repurposing via CMap-style Analysis (fgsea + MSigDB C2:CGP)

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p results/causal_inference/logs
mkdir -p results/drug_repurposing
mkdir -p figures

echo "=== Strategy 11: Drug Repurposing ==="
echo "Start: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-4}"
echo "Memory: $(free -h | head -2)"

# Activate environment
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

# Run the pipeline
Rscript 20_drug_repurposing_cmap.R

echo "End: $(date)"
