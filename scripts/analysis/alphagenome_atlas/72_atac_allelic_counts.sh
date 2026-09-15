#!/usr/bin/env bash
# Step 72 (P3, same-variant chromatin arm): donor-level allele counts at Atlas-scored variants inside snATAC
# consensus peaks, from the GSE281367 cellranger-atac BAMs (hg38, indexed, CB tags present).
# Two on-disk layouts are supported:
#   cellranger : <root>/<donor>/outs/possorted_bam.bam        (GSE281367). cellranger also keeps an internal
#                copy of the same alignment under <donor>/SC_ATAC_COUNTER_CS/.../join-<hash>/files/; matching it
#                too would count every donor twice, so the search is pinned to */outs/.
#   flat       : <root>/<donor>.filtered.bam                  (GSE244832)
# usage: 72_atac_allelic_counts.sh <bam_root> <out_dir> <targets.tsv> [threads] [layout]
set -euo pipefail
BAMROOT=${1:?}; OUT=${2:?}; TARGETS=${3:?}; THREADS=${4:-12}; LAYOUT=${5:-cellranger}
FASTA=/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa
module load BCFtools/1.19-GCC-13.2.0 2>/dev/null || true
mkdir -p "$OUT"
export FASTA TARGETS OUT
one() {
  bam=$1
  if [ "$LAYOUT" = "flat" ]; then donor=$(basename "$bam" | sed 's/\.filtered\.bam$//; s/\.bam$//')
  else donor=$(basename "$(dirname "$(dirname "$bam")")"); fi
  out="$OUT/${donor}.tsv.gz"; [ -s "$out" ] && return 0
  bcftools mpileup -f "$FASTA" -R "$TARGETS" -a AD -d 5000 -q 30 -Q 20 -B \
      --ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY -Ou "$bam" 2>/dev/null \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' 2>/dev/null \
    | awk -v d="$donor" 'BEGIN{OFS="\t"; print "donor","chrom","pos","ref","alt","ad"} {print d,$1,$2,$3,$4,$5}' \
    | gzip -c > "${out}.tmp" && mv "${out}.tmp" "$out"
}
export -f one
export LAYOUT
if [ "$LAYOUT" = "flat" ]; then
  find "$BAMROOT" -maxdepth 1 -name "*.filtered.bam" | sort
else
  find "$BAMROOT" -mindepth 3 -maxdepth 3 -path "*/outs/possorted_bam.bam" | sort
fi | xargs -I{} -P "$THREADS" bash -c 'one "$@"' _ {}
echo "ATAC_ALLELIC_DONE donors=$(ls "$OUT" | wc -l)"
