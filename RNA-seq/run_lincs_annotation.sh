#!/bin/bash
#SBATCH --job-name=liver_lincs_annot
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=results/drug_repurposing/logs/lincs_annotation_%j.out
#SBATCH --error=results/drug_repurposing/logs/lincs_annotation_%j.err

# Script 31: LINCS L1000 Compound Annotation & Integration
# Post-processes LINCS results with compound names, MOA, DGIdb, network proximity
# Run AFTER Script 20 (LINCS query) and optionally after Script 32 (network proximity)

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p results/drug_repurposing/logs
mkdir -p data/lincs_metadata
mkdir -p figures

echo "=== Script 31: LINCS Annotation ==="
echo "Start: $(date)"
echo "Node: $(hostname)"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript 31_lincs_annotation.R

echo "=== Done: $(date) ==="
