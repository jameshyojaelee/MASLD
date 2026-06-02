#!/bin/bash
#SBATCH --job-name=progression_coloc
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/causal_inference/logs/progression_coloc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/causal_inference/logs/progression_coloc_%j.err

# Disease Progression COLOC with Broadaway Liver eQTLs
#
# Usage:
#   GWAS_NAME=GHOUSE_CIRRHOSIS sbatch run_progression_coloc.sh
#   GWAS_NAME=DECODE_NAFL       sbatch run_progression_coloc.sh
#   GWAS_NAME=DECODE_CIRRHOSIS  sbatch run_progression_coloc.sh
#   GWAS_NAME=DECODE_HCC        sbatch run_progression_coloc.sh
#
# Or run Ghouse immediately:
#   GWAS_NAME=GHOUSE_CIRRHOSIS sbatch run_progression_coloc.sh

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
mkdir -p "$BASE/RNA-seq/results/causal_inference/logs"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

GWAS_NAME="${GWAS_NAME:-GHOUSE_CIRRHOSIS}"

echo "============================================================"
echo "  DISEASE PROGRESSION COLOC PIPELINE"
echo "  Started: $(date)"
echo "  SLURM_JOB_ID: $SLURM_JOB_ID"
echo "  GWAS_NAME: $GWAS_NAME"
echo "  RAM: 200G | CPUs: $SLURM_CPUS_PER_TASK"
echo "============================================================"
echo ""

cd "$BASE/RNA-seq"
GWAS_NAME="$GWAS_NAME" Rscript 49_progression_coloc.R 2>&1

echo ""
echo "============================================================"
echo "  COMPLETE: $(date)"
echo "  Results: $BASE/RNA-seq/results/causal_inference/"
echo "============================================================"
