#!/usr/bin/env bash
# Step 30 (P3 feasibility, READ-ONLY inventory): which datasets have per-read alignments or raw reads on disk for
# heterozygous-site calling and allelic counting; which tools exist. Writes p3_feasibility.md + .tsv into $OUT.
set -uo pipefail
P=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
OUT=${OUT:?}
T=$OUT/p3_feasibility.tsv; M=$OUT/p3_feasibility.md
printf "dataset\tpattern\tn_files\ttotal_bytes\texample\n" > $T
inv() {  # dataset, search roots..., pattern
  local ds=$1 pat=$2; shift 2
  local files; files=$(find "$@" -maxdepth 6 -type f \( -iname "$pat" \) 2>/dev/null | head -20000)
  local n; n=$(printf "%s" "$files" | grep -c . || true)
  local bytes=0; if [ "$n" -gt 0 ]; then bytes=$(printf "%s\n" "$files" | head -2000 | xargs -d '\n' stat -c %s 2>/dev/null | awk '{s+=$1} END{print s+0}'); fi
  printf "%s\t%s\t%s\t%s\t%s\n" "$ds" "$pat" "$n" "$bytes" "$(printf "%s\n" "$files" | head -1)" >> $T
}
ROOTS=$(ls -d $P/data $P/RNA-seq $P/Analysis $P/GWAS/finemapping/data /gpfs/commons/groups/sanjana_lab/Cas13 2>/dev/null | tr '\n' ' ')
echo "# P3 feasibility inventory ($(date -u +%FT%TZ))" > $M
for ds in GSE267145 GSE296875 GSE244832 GSE281367 GSE135251 GSE162694 GSE130970 GSE126848 GSE193066 GSE213621 GSE202379; do
  dirs=$(find $ROOTS -maxdepth 5 -type d -iname "*${ds}*" 2>/dev/null | head -10)
  echo "## $ds" >> $M; printf "%s\n" "$dirs" | sed 's/^/- dir: /' >> $M
  if [ -n "$dirs" ]; then
    for pat in "*.bam" "*.cram" "*.fastq.gz" "*.fq.gz" "*.sra" "*fragments.tsv.gz" "possorted*bam" "*.vcf*" "*.bed"; do inv "$ds" "$pat" $dirs; done
  fi
done
echo "## genotype files anywhere" >> $M
find $ROOTS -maxdepth 6 -type f \( -iname "*.vcf" -o -iname "*.vcf.gz" -o -iname "*.bcf" -o -iname "*.bgen" -o -iname "*.pgen" -o -iname "*dosage*" \) 2>/dev/null | grep -v "1kg\|ld_ref\|gsmap" | head -40 | sed 's/^/- /' >> $M
echo "## tools" >> $M
for t in samtools bcftools gatk bedtools STAR bwa; do
  m=$(module avail 2>&1 | grep -io "$t[^ ]*" | head -3 | tr '\n' ' '); e=$(ls -d /gpfs/commons/home/jameslee/micromamba/envs/*/bin/$t 2>/dev/null | head -3 | tr '\n' ' ')
  echo "- $t: modules[$m] envs[$e]" >> $M
done
echo "## reference builds of BAM headers (first BAM per dataset)" >> $M
SAM=$(ls -d /gpfs/commons/home/jameslee/micromamba/envs/*/bin/samtools 2>/dev/null | head -1)
if [ -n "$SAM" ]; then awk -F'\t' 'NR>1 && $2=="*.bam" && $3>0 {print $1"\t"$5}' $T | while IFS=$'\t' read ds f; do echo "- $ds: $($SAM view -H "$f" 2>/dev/null | grep -m1 -o 'SN:chr1[^\t]*\tLN:[0-9]*' ) $($SAM view -H "$f" 2>/dev/null | grep -m1 '@PG' | cut -c1-120)" >> $M; done; fi
echo "DONE" >> $M; echo P3_FEASIBILITY_DONE
