#!/bin/bash
#SBATCH --job-name=liver_net_proximity
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=results/drug_repurposing/logs/network_proximity_%j.out
#SBATCH --error=results/drug_repurposing/logs/network_proximity_%j.err

# Script 32: Network Proximity Scoring (Guney et al. 2016)
# Validates drug repurposing hits via PPI network distance
# Downloads STRING PPI if not present (~500MB)
# 1000 permutations per drug — single-threaded, compute intensive
# If runtime exceeds wall time, reduce N_PERMUTATIONS in 32_network_proximity.R

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p results/drug_repurposing/logs
mkdir -p results/drug_repurposing/network_proximity
mkdir -p data/string_ppi
mkdir -p figures

echo "=== Script 32: Network Proximity Scoring ==="
echo "Start: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-4}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript 32_network_proximity.R

echo "=== Done: $(date) ==="
