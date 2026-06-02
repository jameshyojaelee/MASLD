#!/bin/bash
#SBATCH --job-name=blast_h2m
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --qos=interactive
#SBATCH --output=Cas13_Library_Design/scripts/ortholog_pipeline/blast_work/logs/blast_h2m_%j.out
#SBATCH --error=Cas13_Library_Design/scripts/ortholog_pipeline/blast_work/logs/blast_h2m_%j.err

# Reverse BLAST: human lncRNA query vs mouse lncRNA database.
# Together with the existing mouse->human BLAST this enables reciprocal best
# hit (RBH) filtering -- Wolf & Koonin 2012 (Genome Biol).

set -eo pipefail

PROJECT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
WORK="$PROJECT/Cas13_Library_Design/scripts/ortholog_pipeline/blast_work"
mkdir -p "$WORK/logs"
cd "$PROJECT"

# Activate rnaseq env (provides python) BEFORE module load, then add blast.
# Previous forward job failed because `python` was missing on the PATH at the
# aggregate-step. We need both blastn (via module) and python (via env).
export MAMBA_ROOT_PREFIX="/gpfs/commons/home/jameslee/micromamba"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

module load blast/2.10.0

echo "=== Tooling check ==="
echo "  python : $(command -v python)"
echo "  blastn : $(command -v blastn)"
echo "  makeblastdb : $(command -v makeblastdb)"
python -c "import pandas as pd; print(f'  pandas : {pd.__version__}')"

# Inputs (already extracted by the forward BLAST sbatch)
MOUSE_LNCRNA_FA="$WORK/mouse_vM37_lncRNA.fa"
HUMAN_LNCRNA_FA="$WORK/human_v47_lncRNA.fa"

if [ ! -s "$MOUSE_LNCRNA_FA" ] || [ ! -s "$HUMAN_LNCRNA_FA" ]; then
  echo "ERROR: prerequisite FASTAs missing. Run blast_lncrna_mouse_vs_human.sh first." >&2
  ls -la "$WORK" >&2
  exit 1
fi

echo "  mouse lncRNA FASTA: $(grep -c '^>' $MOUSE_LNCRNA_FA) transcripts"
echo "  human lncRNA FASTA: $(grep -c '^>' $HUMAN_LNCRNA_FA) transcripts"

echo "=== Step 1: Build mouse lncRNA BLAST DB ==="
DB_PREFIX="$WORK/mouse_lncRNA_db"
if [ ! -s "${DB_PREFIX}.nsq" ] && [ ! -s "${DB_PREFIX}.00.nsq" ]; then
  makeblastdb -in "$MOUSE_LNCRNA_FA" -dbtype nucl -out "$DB_PREFIX"
else
  echo "  mouse BLAST DB already exists at $DB_PREFIX.*"
fi

echo "=== Step 2: Run blastn (human lncRNA query vs mouse lncRNA db) ==="
OUT_TSV="$WORK/blast_human_vs_mouse.tsv"
if [ ! -s "$OUT_TSV" ]; then
  blastn \
    -query "$HUMAN_LNCRNA_FA" \
    -db "$DB_PREFIX" \
    -task dc-megablast \
    -evalue 1e-5 \
    -word_size 11 \
    -num_threads ${SLURM_CPUS_PER_TASK:-16} \
    -max_target_seqs 5 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -out "$OUT_TSV"
  echo "  $OUT_TSV :  $(wc -l < $OUT_TSV) lines"
else
  echo "  blast output exists: $(wc -l < $OUT_TSV) lines"
fi

echo "=== Step 3: Aggregate human->mouse transcript hits to gene-level ==="
python "$PROJECT/Cas13_Library_Design/scripts/ortholog_pipeline/aggregate_human_vs_mouse.py" \
  --blast "$OUT_TSV" \
  --mouse-fa "$MOUSE_LNCRNA_FA" \
  --human-fa "$HUMAN_LNCRNA_FA" \
  --out "$WORK/L4_blast_h2m_gene.tsv"

echo "=== Step 4: Enforce reciprocal best hit (RBH) ==="
python "$PROJECT/Cas13_Library_Design/scripts/ortholog_pipeline/enforce_rbh.py" \
  --m2h "$PROJECT/data/external/orthologs/layers/L4_blast.tsv" \
  --h2m "$WORK/L4_blast_h2m_gene.tsv" \
  --out "$PROJECT/data/external/orthologs/layers/L4_blast_rbh.tsv"

echo "=== DONE ==="
date
