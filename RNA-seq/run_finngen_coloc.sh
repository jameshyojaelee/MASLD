#!/bin/bash
#SBATCH --job-name=finngen_coloc
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/causal_inference/logs/finngen_coloc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/causal_inference/logs/finngen_coloc_%j.err

# FinnGen R12 COLOC with Broadaway Liver eQTLs
#
# Usage:
#   GWAS_NAME=FINNGEN_NAFLD sbatch run_finngen_coloc.sh
#   GWAS_NAME=FINNGEN_NASH  sbatch run_finngen_coloc.sh
#   GWAS_NAME=FINNGEN_HCC   sbatch run_finngen_coloc.sh
#
# Or run all three:
#   for g in FINNGEN_NAFLD FINNGEN_NASH FINNGEN_HCC; do
#     GWAS_NAME=$g sbatch run_finngen_coloc.sh
#   done

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
mkdir -p "$BASE/RNA-seq/results/causal_inference/logs"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

GWAS_NAME="${GWAS_NAME:-FINNGEN_NAFLD}"

echo "============================================================"
echo "  FINNGEN R12 COLOC PIPELINE"
echo "  Started: $(date)"
echo "  SLURM_JOB_ID: $SLURM_JOB_ID"
echo "  GWAS_NAME: $GWAS_NAME"
echo "  RAM: 200G | CPUs: $SLURM_CPUS_PER_TASK"
echo "============================================================"
echo ""

# Verify FinnGen data exists
FINNGEN_DIR="$BASE/GWAS/MR_Data/FinnGen"
if [ ! -d "$FINNGEN_DIR" ]; then
    echo "ERROR: FinnGen directory not found: $FINNGEN_DIR"
    exit 1
fi
echo "FinnGen data directory: $(ls -lh "$FINNGEN_DIR"/*.gz 2>/dev/null | wc -l) files"
echo ""

cd "$BASE/RNA-seq"
GWAS_NAME="$GWAS_NAME" Rscript 46_finngen_coloc.R 2>&1

echo ""
echo "============================================================"
echo "  COMPLETE: $(date)"
echo "  Results: $BASE/RNA-seq/results/causal_inference/finngen_*/"
echo "============================================================"
