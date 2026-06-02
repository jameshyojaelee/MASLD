#!/bin/bash
#SBATCH --job-name=bbj_coloc
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/causal_inference/logs/bbj_coloc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/causal_inference/logs/bbj_coloc_%j.err

# BBJ Liver Enzyme COLOC with Broadaway Liver eQTLs (cross-ancestry)
#
# Usage:
#   GWAS_NAME=BBJ_ALT sbatch run_bbj_coloc.sh
#   GWAS_NAME=BBJ_AST sbatch run_bbj_coloc.sh
#   GWAS_NAME=BBJ_GGT sbatch run_bbj_coloc.sh
#
# Or run all three:
#   for g in BBJ_ALT BBJ_AST BBJ_GGT; do
#     GWAS_NAME=$g sbatch run_bbj_coloc.sh
#   done
#
# Prerequisites:
#   1. Download BBJ GWAS: sbatch GWAS/MR_Data/BBJ/download_bbj_liver_enzymes.sh
#   2. Format for COLOC: Rscript GWAS/MR_Data/BBJ/format_bbj_for_coloc.R

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
mkdir -p "$BASE/RNA-seq/results/causal_inference/logs"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

GWAS_NAME="${GWAS_NAME:-BBJ_ALT}"

echo "============================================================"
echo "  BBJ LIVER ENZYME COLOC PIPELINE (cross-ancestry)"
echo "  Started: $(date)"
echo "  SLURM_JOB_ID: $SLURM_JOB_ID"
echo "  GWAS_NAME: $GWAS_NAME"
echo "  GWAS ancestry: East Asian (Japanese)"
echo "  eQTL ancestry: European (Broadaway N=1,183)"
echo "  RAM: 200G | CPUs: $SLURM_CPUS_PER_TASK"
echo "============================================================"
echo ""

# Verify BBJ harmonised data exists
BBJ_FILE="$BASE/GWAS/MR_Data/BBJ/${GWAS_NAME#BBJ_}"
BBJ_HARMONISED="$BASE/GWAS/MR_Data/BBJ/BBJ_${GWAS_NAME#BBJ_}_harmonised_hg38.tsv.gz"
if [ ! -f "$BBJ_HARMONISED" ]; then
    echo "ERROR: Harmonised BBJ file not found: $BBJ_HARMONISED"
    echo "Run format_bbj_for_coloc.R first:"
    echo "  Rscript GWAS/MR_Data/BBJ/format_bbj_for_coloc.R"
    exit 1
fi
echo "BBJ harmonised file: $(ls -lh "$BBJ_HARMONISED" | awk '{print $5}')"
echo ""

cd "$BASE/RNA-seq"
GWAS_NAME="$GWAS_NAME" Rscript 47_bbj_coloc.R 2>&1

echo ""
echo "============================================================"
echo "  COMPLETE: $(date)"
echo "  Results: $BASE/RNA-seq/results/causal_inference/bbj_*/"
echo "============================================================"
