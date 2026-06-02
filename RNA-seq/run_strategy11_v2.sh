#!/bin/bash
#SBATCH --job-name=liver_drug_repurpose_v2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=results/causal_inference/logs/drug_repurposing_v2_%j.out
#SBATCH --error=results/causal_inference/logs/drug_repurposing_v2_%j.err

# Strategy 11 v2: Drug Repurposing with t-stat ranking + signatureSearch (if available)
# Key improvements:
#   - t-statistic ranking (not logFC)
#   - signatureSearch LINCS L1000 query (HepG2) if installed
#   - Increased nPermSimple=10000 for stable p-values
#   - CGP results reframed as disease concordance

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p results/causal_inference/logs
mkdir -p results/drug_repurposing
mkdir -p figures

echo "=== Strategy 11 v2: Drug Repurposing ==="
echo "Start: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-4}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript 20_drug_repurposing_cmap.R

echo "=== Done: $(date) ==="
