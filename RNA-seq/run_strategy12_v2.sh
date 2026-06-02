#!/bin/bash
#SBATCH --job-name=liver_gwas_overlay_v2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=results/gwas_spatial_convergence/logs/gwas_overlay_v2_%j.out
#SBATCH --error=results/gwas_spatial_convergence/logs/gwas_overlay_v2_%j.err

# Strategy 12 v2: GWAS Overlay with Human-Only Mapping + Continuous Enrichment

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p results/gwas_spatial_convergence/logs
mkdir -p figures

echo "=== Strategy 12 v2: GWAS Overlay ==="
echo "Start: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-4}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

Rscript 30_gwas_overlay_v2.R

echo "=== Done: $(date) ==="
