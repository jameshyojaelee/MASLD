#!/usr/bin/env bash
#SBATCH --job-name=alphagenome
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/79_alphagenome_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/79_alphagenome_%j.err

# Phase-5 src/79 AlphaGenome panel: caQTL ACCESSIBILITY (ATAC+DNASE heads) + sQTL
# SPLICE (SpliceJunction+SpliceSiteUsage) over the leads that lack features. No GPU
# (Google-side); only gRPC egress. Key from $ALPHA_GENOME_API_KEY / ~/.alphagenome_key,
# never logged. Modality-matched. ADDITIVE; edits nothing.
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRC="${ROOT}/GWAS/finemapping/src"
SEQFUNC="${ROOT}/GWAS/finemapping/results/seqfunc"
mkdir -p "${SEQFUNC}/logs"

echo "=== src/79 AlphaGenome caQTL-accessibility + sQTL-splice ==="; echo "$(hostname) | $(date)"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate alphagenome
export PYTHONUNBUFFERED=1 PYTHONNOUSERSITE=1

WHICH="${WHICH:-both}"
SEQLEN="${SEQLEN:-SEQUENCE_LENGTH_500KB}"
MAXVARS="${MAXVARS:-0}"      # API can handle thousands; 0=all
python3 "${SRC}/79_score_panels.py" --panel alphagenome \
    --which "${WHICH}" --seq-length "${SEQLEN}" --max-variants "${MAXVARS}"
echo "Done: $(date)"
