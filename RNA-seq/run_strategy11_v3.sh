#!/bin/bash
#SBATCH --job-name=liver_drug_repurpose_v3
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=results/drug_repurposing/logs/strategy11_v3_%j.out
#SBATCH --error=results/drug_repurposing/logs/strategy11_v3_%j.err

# Strategy 11 v3: Comprehensive Pharmacotranscriptomics (with critical fixes)
#   - Arm 1: LINCS L1000 signature reversal (HepG2, signatureSearch)
#     + BRD→drug name mapping via CLUE Repurposing Hub (downloads compoundinfo_beta.txt)
#     + WTCS_FDR significance filter replacing Tau
#   - Arm 2a: DGIdb drug-gene interactions (GraphQL API)
#   - Arm 2b: Open Targets known drugs (expanded: MR-causal + DGIdb-druggable genes)
#   - C2:CGP disease concordance (fgsea) with concordant leading-edge for upregulated genes
#   - Integration + MR convergence + composite figure
# Requires internet: CLUE S3, DGIdb API, Open Targets API

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p results/drug_repurposing/logs
mkdir -p figures

echo "=== Strategy 11 v3: Comprehensive Pharmacotranscriptomics ==="
echo "Start: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-8}"
echo "Memory: ${SLURM_MEM_PER_NODE:-64G}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

export SLURM_CPUS=${SLURM_CPUS_PER_TASK:-8}

Rscript 20_drug_repurposing_v3.R

echo "=== Done: $(date) ==="
