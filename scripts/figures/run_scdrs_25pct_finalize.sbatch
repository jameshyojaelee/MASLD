#!/bin/bash
#SBATCH --job-name=scdrs_25pct_fin
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=2:00:00
#SBATCH --output=scripts/figures/logs/scdrs_25pct_finalize_%j.out
#SBATCH --error=scripts/figures/logs/scdrs_25pct_finalize_%j.err

# Aggregate the 25% subsample variant + regenerate figS_scdrs_gwas.pdf.

set -eo pipefail
mkdir -p scripts/figures/logs
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"

echo "[finalize 25pct] step 1: aggregate"
micromamba activate spatial
python scripts/figures/figS_scdrs_null_aggregate_25pct.py
micromamba deactivate

echo "[finalize 25pct] step 2: regenerate figure"
micromamba activate rnaseq
Rscript scripts/figures/figS_scdrs_gwas.R

echo "[finalize 25pct] done. Hepatocyte summary:"
awk -F, 'NR==1 || $2=="Hepatocytes"' Analysis/SingleCell/results_gpu_v2/disease_signatures/scdrs_gwas_hep_null_summary_25pct.csv
