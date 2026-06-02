#!/bin/bash
#SBATCH --job-name=blast_lncrna
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=Cas13_Library_Design/scripts/ortholog_pipeline/blast_work/logs/blast_%j.out
#SBATCH --error=Cas13_Library_Design/scripts/ortholog_pipeline/blast_work/logs/blast_%j.err

# BLAST mouse lncRNA transcripts vs human lncRNA database.
# Output is L4 (sequence-level) layer for the ortholog table.

set -eo pipefail

PROJECT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
WORK="$PROJECT/Cas13_Library_Design/scripts/ortholog_pipeline/blast_work"
mkdir -p "$WORK/logs"
cd "$PROJECT"

# Activate rnaseq env (provides python + pandas) BEFORE module load.
# Previous run failed at Step 4 because python was not on PATH.
export MAMBA_ROOT_PREFIX="/gpfs/commons/home/jameslee/micromamba"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

module load blast/2.10.0

echo "=== Tooling check ==="
echo "  python : $(command -v python)"
echo "  blastn : $(command -v blastn)"
python -c "import pandas as pd; print(f'  pandas : {pd.__version__}')"

# Inputs
HUMAN_FA="/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/genome/gencode.v47.lncRNA_transcripts.fa.gz"
MOUSE_FA="/gpfs/commons/groups/sanjana_lab/Cas13/TIGER/references/gencode.vM37.transcripts.fa"

echo "=== Step 1: Extract mouse lncRNA-only FASTA with simplified headers ==="
MOUSE_LNCRNA_FA="$WORK/mouse_vM37_lncRNA.fa"
MOUSE_LNCRNA_HEADERS="$WORK/mouse_vM37_lncRNA_headers.tsv"
if [ ! -s "$MOUSE_LNCRNA_FA" ] || [ ! -s "$MOUSE_LNCRNA_HEADERS" ]; then
  # GENCODE header: >tx_id|gene_id|otg_id|ottt_id|tx_name|gene_name|length|biotype|
  # Keep biotype == lncRNA. Simplify header to just tx_id (no pipes) and dump
  # full mapping to a side TSV.
  rm -f "$MOUSE_LNCRNA_FA" "$MOUSE_LNCRNA_HEADERS"
  echo -e "tx_id\tgene_id\tgene_name\tlength" > "$MOUSE_LNCRNA_HEADERS"
  awk -v hdr="$MOUSE_LNCRNA_HEADERS" 'BEGIN{keep=0}
       /^>/ {
         keep=0
         line=substr($0,2)
         n=split(line, f, "|")
         if (n>=8 && f[8]=="lncRNA") {
           keep=1
           print ">" f[1]
           print f[1] "\t" f[2] "\t" f[6] "\t" f[7] >> hdr
         }
       }
       !/^>/ && keep {print}' "$MOUSE_FA" > "$MOUSE_LNCRNA_FA"
  echo "  $MOUSE_LNCRNA_FA :  $(grep -c '^>' $MOUSE_LNCRNA_FA) transcripts"
fi

echo "=== Step 2: Build human lncRNA FASTA with simplified headers + BLAST DB ==="
HUMAN_LNCRNA_FA="$WORK/human_v47_lncRNA.fa"
HUMAN_LNCRNA_HEADERS="$WORK/human_v47_lncRNA_headers.tsv"
if [ ! -s "$HUMAN_LNCRNA_FA" ] || [ ! -s "$HUMAN_LNCRNA_HEADERS" ]; then
  rm -f "$HUMAN_LNCRNA_FA" "$HUMAN_LNCRNA_HEADERS"
  echo -e "tx_id\tgene_id\tgene_name\tlength" > "$HUMAN_LNCRNA_HEADERS"
  zcat "$HUMAN_FA" | awk -v hdr="$HUMAN_LNCRNA_HEADERS" '/^>/ {
         line=substr($0,2)
         n=split(line, f, "|")
         print ">" f[1]
         print f[1] "\t" f[2] "\t" f[6] "\t" f[7] >> hdr
         next
       }
       {print}' > "$HUMAN_LNCRNA_FA"
  echo "  $HUMAN_LNCRNA_FA :  $(grep -c '^>' $HUMAN_LNCRNA_FA) transcripts"
fi

DB_PREFIX="$WORK/human_lncRNA_db"
if [ ! -s "${DB_PREFIX}.nsq" ] && [ ! -s "${DB_PREFIX}.00.nsq" ]; then
  makeblastdb -in "$HUMAN_LNCRNA_FA" -dbtype nucl -out "$DB_PREFIX"
fi

echo "=== Step 3: Run blastn (mouse lncRNA query vs human lncRNA db) ==="
OUT_TSV="$WORK/blast_mouse_vs_human.tsv"
if [ ! -s "$OUT_TSV" ]; then
  blastn \
    -query "$MOUSE_LNCRNA_FA" \
    -db "$DB_PREFIX" \
    -task dc-megablast \
    -evalue 1e-5 \
    -word_size 11 \
    -num_threads ${SLURM_CPUS_PER_TASK:-16} \
    -max_target_seqs 50 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -out "$OUT_TSV"
  echo "  $OUT_TSV :  $(wc -l < $OUT_TSV) lines"
else
  echo "  blast output exists: $(wc -l < $OUT_TSV) lines"
fi

echo "=== Step 4: Aggregate transcript-level hits to gene-level pairs ==="
python "$PROJECT/Cas13_Library_Design/scripts/ortholog_pipeline/blast_to_orthotable.py" \
  --blast "$OUT_TSV" \
  --mouse-fa "$MOUSE_LNCRNA_FA" \
  --human-fa "$HUMAN_LNCRNA_FA" \
  --out "$PROJECT/data/external/orthologs/layers/L4_blast.tsv"

echo "=== DONE ==="
date
