#!/usr/bin/env bash
# Download and prepare GRCm39 references and indices.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
REF_DIR="$ROOT/reference"
RAW_DIR="$REF_DIR/raw"
STAR_DIR="$REF_DIR/indices/star"
STAR2_DIR="$REF_DIR/indices/star_2pass"
HISAT2_DIR="$REF_DIR/indices/hisat2"
SALMON_DIR="$REF_DIR/indices/salmon"
THREADS=${THREADS:-16}
GENCODE_BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33"
FASTA_NAME="GRCm39.primary_assembly.genome.fa"
GTF_NAME="gencode.vM33.annotation.gtf"
TRANSCRIPTS_FA="gencode.vM33.transcripts.fa"
RSEQC_BED="gencode.vM33.annotation.rseqc.bed"
mkdir -p "$RAW_DIR" "$STAR_DIR" "$STAR2_DIR" "$HISAT2_DIR" "$SALMON_DIR"
cd "$RAW_DIR"
if [[ ! -s "$FASTA_NAME" ]]; then
  wget -O "${FASTA_NAME}.gz" "$GENCODE_BASE/${FASTA_NAME}.gz"
  gunzip -f "${FASTA_NAME}.gz"
fi
if [[ ! -s "$GTF_NAME" ]]; then
  wget -O "${GTF_NAME}.gz" "$GENCODE_BASE/${GTF_NAME}.gz"
  gunzip -f "${GTF_NAME}.gz"
fi
if [[ ! -s "$TRANSCRIPTS_FA" ]]; then
  gffread "$GTF_NAME" -g "$FASTA_NAME" -w "$TRANSCRIPTS_FA"
fi
if [[ ! -s "$RSEQC_BED" ]]; then
  python - "$RAW_DIR" "$GTF_NAME" "$RSEQC_BED" <<'PY'
import sys
from pathlib import Path
from collections import defaultdict
raw_dir = Path(sys.argv[1])
gtf_path = raw_dir / sys.argv[2]
bed_path = raw_dir / sys.argv[3]
features = defaultdict(list)
with gtf_path.open() as handle:
    for line in handle:
        if line.startswith('#'):
            continue
        chrom, source, feature, start, end, score, strand, frame, attrs = line.strip().split('\t')
        if feature != 'exon':
            continue
        attr_dict = {}
        for field in attrs.split(';'):
            field = field.strip()
            if not field:
                continue
            key, value = field.split(' ', 1)
            attr_dict[key] = value.strip('"')
        transcript_id = attr_dict.get('transcript_id')
        gene_id = attr_dict.get('gene_id', transcript_id)
        if not transcript_id:
            continue
        features[transcript_id].append((chrom, int(start) - 1, int(end), strand, gene_id))
records = []
for transcript_id, exons in features.items():
    exons.sort(key=lambda x: x[1])
    chrom = exons[0][0]
    strand = exons[0][3]
    gene_id = exons[0][4]
    starts = [exon[1] for exon in exons]
    ends = [exon[2] for exon in exons]
    tx_start = min(starts)
    tx_end = max(ends)
    block_sizes = [str(end - start) for start, end in zip(starts, ends)]
    block_starts = [str(start - tx_start) for start in starts]
    name = f"{gene_id}|{transcript_id}"
    record = (chrom, tx_start, tx_end, name, '0', strand, tx_start, tx_start,
              '0', len(exons), ','.join(block_sizes), ','.join(block_starts))
    records.append('\t'.join(map(str, record)))
bed_path.write_text('\n'.join(records) + '\n')
PY
fi
cd "$STAR_DIR"
if [[ ! -f Genome ]]; then
  STAR --runThreadN "$THREADS" \
       --runMode genomeGenerate \
       --genomeFastaFiles "$RAW_DIR/$FASTA_NAME" \
       --sjdbGTFfile "$RAW_DIR/$GTF_NAME" \
       --sjdbOverhang 149
fi
cd "$SALMON_DIR"
if [[ ! -f versionInfo.json ]]; then
  salmon index -t "$RAW_DIR/$TRANSCRIPTS_FA" -i "$SALMON_DIR" -k 31
fi
cd "$HISAT2_DIR"
if [[ ! -f grcm39.1.ht2 ]]; then
  hisat2-build "$RAW_DIR/$FASTA_NAME" "$HISAT2_DIR/grcm39"
fi
echo "Reference setup complete."
