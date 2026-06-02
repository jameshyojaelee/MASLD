#!/bin/bash
#SBATCH --job-name=magma_fin
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=2:00:00
#SBATCH --output=scripts/figures/logs/magma_fin_%j.out
#SBATCH --error=scripts/figures/logs/magma_fin_%j.err

# Runs after the scdrs_magma scoring array completes (--dependency=afterok).
# Steps: (a) magma_to_gs.py — convert MAGMA .genes.out → scDRS .gs (already
# done before scoring, but harmless to run again as no-op since .gs exists);
# (b) figS_scdrs_magma_aggregate.py — combine per-GWAS CT enrichment;
# (c) figS_scdrs_magma.R — build figS_scdrs_magma.pdf with 4 panels.

set -eo pipefail
mkdir -p scripts/figures/logs
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"

echo "[finalize] aggregate"
micromamba activate spatial
python scripts/figures/figS_scdrs_magma_aggregate.py
micromamba deactivate

echo "[finalize] render figure"
micromamba activate rnaseq
Rscript scripts/figures/figS_scdrs_magma.R

echo "[finalize] done. Top-CT summary per GWAS:"
awk -F, 'NR<=20' Analysis/SingleCell/results_gpu_v2/disease_signatures/scdrs_magma_summary.csv 2>/dev/null
