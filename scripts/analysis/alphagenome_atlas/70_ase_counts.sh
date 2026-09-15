#!/usr/bin/env bash
# Step 70 (P3, bulk-RNA arm): allele counts at Resource variants from the cohort STAR BAMs.
#
# The BAMs are hg38 with chr naming and are NOT indexed, so bcftools streams each file once with -T (targets),
# which needs no index. Only the prespecified target sites are counted; no genotypes are called here.
# usage: 70_ase_counts.sh <cohort> <bam_dir> <out_dir> <targets.tsv.gz> [threads]
set -euo pipefail
COHORT=${1:?}; BAMDIR=${2:?}; OUT=${3:?}; TARGETS=${4:?}; THREADS=${5:-16}
FASTA=/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa
module load BCFtools/1.19-GCC-13.2.0 2>/dev/null || true
command -v bcftools >/dev/null || { echo "bcftools not on PATH"; exit 1; }
mkdir -p "$OUT/counts"
export FASTA TARGETS OUT
count_one() {
  bam=$1
  sample=$(basename "$bam" | sed 's/\.Aligned.*//; s/_Aligned.*//; s/\.bam$//')
  out="$OUT/counts/${sample}.tsv.gz"
  [ -s "$out" ] && return 0
  # -q 10 MAPQ, -Q 20 base quality, -B no BAQ (RNA), --ff drop dup/secondary/supplementary/unmapped/QC-fail
  bcftools mpileup -f "$FASTA" -T "$TARGETS" -a AD -d 1000 -q 10 -Q 20 -B \
      --ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY -Ou "$bam" 2>/dev/null \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' 2>/dev/null \
    | awk -v s="$sample" 'BEGIN{OFS="\t"; print "sample","chrom","pos","ref","alt","ad"} {print s,$1,$2,$3,$4,$5}' \
    | gzip -c > "${out}.tmp" && mv "${out}.tmp" "$out"
}
export -f count_one
find "$BAMDIR" -name "*.bam" | sort | xargs -I{} -P "$THREADS" bash -c 'count_one "$@"' _ {}
n=$(ls "$OUT/counts" | wc -l)
echo "ASE_COUNTS_DONE $COHORT samples=$n"
