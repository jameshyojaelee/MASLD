#!/usr/bin/env bash
# D1 task 2: per donor x lineage allele counts at the 1,119 P3 target sites.
#
# The pileup recipe is copied unchanged from scripts/analysis/alphagenome_atlas/72_atac_allelic_counts.sh
# (bcftools mpileup -a AD -d 5000 -q 30 -Q 20 -B --ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY, the
# cellranger-arc GRCh38-2024-A FASTA). That script is not modified. What is added here is the lineage
# split: the donor BAM is first reduced to reads overlapping the 1,119 targets (samtools -M -L, index
# driven), then one stream per lineage is selected by barcode, then each stream is piled up separately.
#
# Two barcode conventions:
#   cb   GSE281367, cellranger-atac: the CB:Z tag.
#   name GSE244832, combinatorial indexing: the read-name prefix before the first ':'.
#
# ALL_READS is the no-barcode-filter stream. It reproduces 72_atac_allelic_counts.sh and exists only so
# the deposit's site count can be regenerated through this loader.
#
# usage: d1_02_extract_and_count.sh <out_root> <threads>
set -euo pipefail
OUT=${1:?out_root}; THREADS=${2:-16}
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
FASTA=/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa
TARGETS=$PROJECT/GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/tables/atac_targets_in_peaks.tsv
BCDIR=$OUT/barcodes
SUBDIR=$OUT/raw/subset_bam
CNTDIR=$OUT/raw/counts
QCDIR=$OUT/raw/qc
mkdir -p "$SUBDIR" "$CNTDIR" "$QCDIR"

BED=$OUT/raw/targets.bed
awk 'BEGIN{OFS="\t"} {print $1,$2-1,$2}' "$TARGETS" | sort -k1,1 -k2,2n > "$BED"

LINEAGES="Hepatocyte Endothelial Stellate_Cell Macrophage Kupffer_Cell Cholangiocyte LSEC Plasma_Cell NK_T_Cell B_Cell Myeloid ALL_LABELLED"

export OUT PROJECT FASTA TARGETS BCDIR SUBDIR CNTDIR QCDIR BED LINEAGES

one_donor() {
  set -euo pipefail
  cohort=$1; donor=$2; bam=$3; mode=$4
  sub="$SUBDIR/${donor}.bam"
  if [ ! -s "$sub" ]; then
    samtools view -b -M -L "$BED" -o "${sub}.tmp.bam" "$bam"
    samtools index "${sub}.tmp.bam"
    mv "${sub}.tmp.bam" "$sub"; mv "${sub}.tmp.bam.bai" "${sub}.bai"
  fi
  total=$(samtools view -c "$sub")
  echo -e "${cohort}\t${donor}\tALL_READS\t${total}" > "$QCDIR/${donor}.reads.tsv"

  pileup() {  # pileup <bam> <lineage>
    b=$1; lin=$2
    out="$CNTDIR/${donor}__${lin}.tsv.gz"
    [ -s "$out" ] && return 0
    # a lineage stream with no read at any target is a real result, not an error: emit a header-only
    # file so the assembler sees the stream and one failing stream never aborts the other 29 donors
    if ! { bcftools mpileup -f "$FASTA" -R "$TARGETS" -a AD -d 5000 -q 30 -Q 20 -B \
             --ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY -Ou "$b" 2>/dev/null \
           | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' 2>/dev/null \
           | awk -v d="$donor" -v l="$lin" -v c="$cohort" 'BEGIN{OFS="\t"; print "cohort","donor","lineage","chrom","pos","ref","alt","ad"} {print c,d,l,$1,$2,$3,$4,$5}' \
           | gzip -c > "${out}.tmp"; }; then
      echo "D1_PILEUP_FAIL ${donor} ${lin}" >&2
      printf 'cohort\tdonor\tlineage\tchrom\tpos\tref\talt\tad\n' | gzip -c > "${out}.tmp"
    fi
    mv "${out}.tmp" "$out"
  }

  pileup "$sub" ALL_READS

  for lin in $LINEAGES; do
    bcl="$BCDIR/${donor}__${lin}.txt"
    [ -s "$bcl" ] || { echo -e "${cohort}\t${donor}\t${lin}\t0" >> "$QCDIR/${donor}.reads.tsv"; continue; }
    lb="$SUBDIR/${donor}__${lin}.bam"
    samtools view -h "$sub" \
      | awk -v F="$bcl" -v MODE="$mode" '
          BEGIN{while((getline l < F) > 0) keep[l]=1}
          /^@/{print; next}
          { bc="";
            if (MODE=="cb") { for(i=12;i<=NF;i++) if (substr($i,1,5)=="CB:Z:") { bc=substr($i,6); break } }
            else { n=split($1,a,":"); bc=a[1] }
            if (bc in keep) print }' \
      | samtools view -b -o "$lb" -
    samtools index "$lb"
    n=$(samtools view -c "$lb")
    echo -e "${cohort}\t${donor}\t${lin}\t${n}" >> "$QCDIR/${donor}.reads.tsv"
    if [ "$n" -gt 0 ]; then pileup "$lb" "$lin"; else
      printf 'cohort\tdonor\tlineage\tchrom\tpos\tref\talt\tad\n' | gzip -c > "$CNTDIR/${donor}__${lin}.tsv.gz"
    fi
    rm -f "$lb" "${lb}.bai"
  done
  echo "D1_DONOR_DONE ${donor} reads=${total}"
}
export -f one_donor

MANIFEST=$OUT/raw/donor_manifest.tsv
: > "$MANIFEST"
for d in $(seq -w 1 12); do
  b="$PROJECT/Analysis/ATAC/Human_External/cellranger/Z${d}/outs/possorted_bam.bam"
  [ -s "$b" ] && echo -e "GSE281367\tZ${d}\t${b}\tcb" >> "$MANIFEST"
done
for d in $(seq -w 1 18); do
  b="$PROJECT/Analysis/ATAC/Human_Multiome/results/alignment/D${d}.filtered.bam"
  [ -s "$b" ] && echo -e "GSE244832\tD${d}\t${b}\tname" >> "$MANIFEST"
done
echo "D1 donors found: $(wc -l < "$MANIFEST")"

# no path in the manifest contains a space, so a space-separated xargs -n 4 dispatch is safe and readable
tr '\t' ' ' < "$MANIFEST" | xargs -n 4 -P "$THREADS" bash -c 'one_donor "$0" "$1" "$2" "$3"'

cat "$QCDIR"/*.reads.tsv | awk 'BEGIN{OFS="\t"; print "cohort","donor","lineage","n_reads_at_targets"} {print}' > "$OUT/tables/reads_at_targets.tsv"
echo "D1_COUNTS_DONE files=$(ls "$CNTDIR" | wc -l)"
