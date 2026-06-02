#!/bin/bash
#SBATCH --job-name=scdrs_null_finalize
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=4:00:00
#SBATCH --output=scripts/figures/logs/scdrs_null_finalize_%j.out
#SBATCH --error=scripts/figures/logs/scdrs_null_finalize_%j.err

# Final stage of the bias-controlled scDRS null pipeline. Runs after the
# 200-task SLURM null-scoring array completes via --dependency=afterok.
# Steps:
#   1. Aggregate per-null per-CT enrichment CSVs → empirical-p + Δz summary
#   2. Re-render figS_scdrs_gwas.pdf with Panel A2 (bias-controlled null)
#   3. Re-render the standalone panel PDFs (figS_scdrs_gwas_A2_*.pdf will appear)

set -eo pipefail
mkdir -p scripts/figures/logs

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"

echo "[finalize] step 1: aggregate null + real CT enrichment"
micromamba activate spatial
python scripts/figures/figS_scdrs_null_aggregate.py
micromamba deactivate

echo "[finalize] step 2: regenerate figS_scdrs_gwas.pdf with Panel A2"
micromamba activate rnaseq
Rscript scripts/figures/figS_scdrs_gwas.R

echo "[finalize] done. Outputs:"
ls -la figures/supplementary/figS_scdrs/figS_scdrs_gwas.pdf
ls figures/supplementary/figS_scdrs/panels/ | grep A2
