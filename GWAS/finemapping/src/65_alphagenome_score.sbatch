#!/usr/bin/env bash
#SBATCH --job-name=alphagenome
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/alphagenome_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/alphagenome_%j.err

# ---------------------------------------------------------------------------
# Stage-2b seqfunc: AlphaGenome API scorer over the 234 independent Broadaway
# eQTL lead signals (+ SORT1 / PNPLA3 / HSD17B13 anchors). Independent SECOND
# sequence-to-function model alongside Borzoi. Additive-only (script 65).
#
# NO GPU: AlphaGenome inference is Google-side; this job only opens a gRPC
# channel to gdmscience.googleapis.com:443. Runs in the persistent `alphagenome`
# micromamba env (python 3.11 + alphagenome 0.7.0 + pysam).
#
# API key: read from $ALPHA_GENOME_API_KEY, else ~/.alphagenome_key (chmod 600).
# The key is never logged and never written into the repo. If no key is present,
# the script prints a READY-awaiting-key message and exits 3.
# ---------------------------------------------------------------------------
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SEQFUNC="${ROOT}/GWAS/finemapping/results/seqfunc"
SRC="${ROOT}/GWAS/finemapping/src"

FASTA=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa
TRUTH="${SEQFUNC}/broadaway_benchmark_truth.tsv"
BORZOI="${SEQFUNC}/borzoi_eqtl_scores.tsv"
OUT="${SEQFUNC}/alphagenome_eqtl_scores.tsv"

mkdir -p "${SEQFUNC}/logs"

echo "=== AlphaGenome API eQTL scorer ==="
echo "Node: $(hostname) | Date: $(date)"
echo ""

# --- env (isolated from ~/.local user-site; the alphagenome env leaks otherwise) ---
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate alphagenome
export PYTHONUNBUFFERED=1
export PYTHONNOUSERSITE=1

python3 "${SRC}/65_alphagenome_score.py" \
    --truth         "${TRUTH}" \
    --fasta         "${FASTA}" \
    --borzoi-scores "${BORZOI}" \
    --seq-length    SEQUENCE_LENGTH_1MB \
    --out           "${OUT}"

echo "Done: $(date)"
