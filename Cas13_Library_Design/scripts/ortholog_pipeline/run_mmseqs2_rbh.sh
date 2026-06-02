#!/usr/bin/env bash
#SBATCH --job-name=mmseqs2_rbh
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=mmseqs2_work/logs/mmseqs2_rbh_%j.out
#SBATCH --error=mmseqs2_work/logs/mmseqs2_rbh_%j.err
set -euo pipefail

# ── MMseqs2 easy-rbh for lncRNA ortholog mapping (L5) ───────────────────────
# Runs bidirectional nucleotide search + reciprocal best hit in one command.
# Input : pre-simplified lncRNA FASTAs (from BLAST pipeline)
# Output: BLAST-like tab file → parsed by parse_mmseqs2.py into L5 layer

module purge
module load bio/MMseqs2/18-8cc5c-gompi-2025a

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
BLAST_WORK="${PIPELINE_DIR}/blast_work"
WORK="${PIPELINE_DIR}/mmseqs2_work"
mkdir -p "${WORK}/tmp" "${WORK}/logs"

MOUSE_FA="${BLAST_WORK}/mouse_vM37_lncRNA.fa"
HUMAN_FA="${BLAST_WORK}/human_v47_lncRNA.fa"

echo "[$(date)] MMseqs2 version: $(mmseqs version 2>&1 || echo 'unknown')"
echo "[$(date)] Mouse FASTA: $(grep -c '^>' "${MOUSE_FA}") sequences"
echo "[$(date)] Human FASTA: $(grep -c '^>' "${HUMAN_FA}") sequences"
echo "[$(date)] Threads: ${SLURM_CPUS_PER_TASK:-16}"

# ── Step 1: Run easy-rbh (nucleotide mode) ──────────────────────────────────
# --search-type 3 = nucleotide-vs-nucleotide
# -e 1e-5         = E-value threshold (matching BLAST pipeline)
# --min-seq-id 0.5 = minimum sequence identity 50% (matching BLAST pipeline)
# -s 1             = fast prefilter (lowest practical sensitivity; the 152K x
#                    191K nucl prefilter at default k=6, exact-kmer-matching is
#                    extremely slow at any s; s=1 minimizes the diagonal scoring
#                    threshold). Combined with --mask 1 to skip low-complexity
#                    homopolymer repeats that dominate the k-mer index.
# --exact-kmer-matching 0 = allow approximate k-mer matching (huge speedup for
#                    nucleotide mode where default is exact-only)
# Output columns match BLAST outfmt 6

echo "[$(date)] Starting MMseqs2 easy-rbh..."

# Clean previous tmp to avoid stale DB conflicts
rm -rf "${WORK}/tmp"/*

mmseqs easy-rbh \
    "${MOUSE_FA}" \
    "${HUMAN_FA}" \
    "${WORK}/mmseqs2_rbh_result" \
    "${WORK}/tmp" \
    --threads "${SLURM_CPUS_PER_TASK:-16}" \
    --search-type 3 \
    -s 1 \
    -e 1e-5 \
    --min-seq-id 0.5 \
    --mask 1 \
    --exact-kmer-matching 0 \
    --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

echo "[$(date)] MMseqs2 easy-rbh complete."

# ── Step 2: Basic stats ─────────────────────────────────────────────────────
RESULT="${WORK}/mmseqs2_rbh_result"
if [[ -f "${RESULT}" ]]; then
    echo "[$(date)] Output lines: $(wc -l < "${RESULT}")"
    echo "[$(date)] Unique query IDs: $(cut -f1 "${RESULT}" | sort -u | wc -l)"
    echo "[$(date)] Unique target IDs: $(cut -f2 "${RESULT}" | sort -u | wc -l)"
    echo "[$(date)] First 5 lines:"
    head -5 "${RESULT}"
else
    echo "[$(date)] ERROR: output file not found at ${RESULT}"
    exit 1
fi

echo "[$(date)] Done. Next: python parse_mmseqs2.py"
