#!/bin/bash
#SBATCH --job-name=panukbb_coloc
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/causal_inference/logs/panukbb_coloc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/causal_inference/logs/panukbb_coloc_%j.err

# Pan-UKBB Multi-Ancestry COLOC
#
# Usage:
#   GWAS_NAME=PANUKBB_AFR_ALT sbatch run_panukbb_coloc.sh
#   GWAS_NAME=PANUKBB_CSA_GGT sbatch run_panukbb_coloc.sh
#
# Run all 6 (AFR + CSA x ALT/AST/GGT):
#   for pop in AFR CSA; do for trait in ALT AST GGT; do
#     GWAS_NAME=PANUKBB_${pop}_${trait} sbatch run_panukbb_coloc.sh
#   done; done

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
mkdir -p "$BASE/RNA-seq/results/causal_inference/logs"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

GWAS_NAME="${GWAS_NAME:-PANUKBB_AFR_ALT}"

echo "============================================================"
echo "  PAN-UKBB MULTI-ANCESTRY COLOC"
echo "  Started: $(date)"
echo "  SLURM_JOB_ID: $SLURM_JOB_ID"
echo "  GWAS_NAME: $GWAS_NAME"
echo "============================================================"

cd "$BASE/RNA-seq"
GWAS_NAME="$GWAS_NAME" Rscript 50_panukbb_coloc.R 2>&1

echo "COMPLETE: $(date)"
