#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --array=0-219
# A4 step 1: RNA genotype calls, one cohort x one autosome per task (10 x 22),
# on the existing STAR BAMs through the A1 symlinked indexes, over merged GENCODE
# v49 exons (PRESPEC_A section 6): MAPQ 255, BQ 20, depth cap 1,000, SNVs,
# biallelic, QUAL >= 30, REDIportal / ENCODE blacklist / MHC removed, DP < 10 -> missing.
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%A_%a.out a4_01_call.sbatch <a1_restricted_dir> <a4_restricted_dir>
set -euo pipefail
a1=$1; out=$2
cohorts=(GSE126848 GSE130970 GSE135251 GSE162694 GSE167523 GSE174478 GSE193066 GSE213621 GSE240729 PRJNA512027)
cohort=${cohorts[$((SLURM_ARRAY_TASK_ID / 22))]}
c=chr$((SLURM_ARRAY_TASK_ID % 22 + 1))
ref=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa
gtf=/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49/gencode.v49.primary_chr.gtf
ext=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external/allelic_refs
bl=$(ls /gpfs/commons/home/jameslee/reference_genome/blacklists/*hg38*.bed* 2>/dev/null | awk 'NR == 1')
module load bcftools/1.21 htslib/1.21 >/dev/null 2>&1
mkdir -p "$out/calls" "$out/regions"

reg=$out/regions/$c.exons_minus_excluded.bed
if [ ! -s "$reg" ]; then
  awk -F'\t' -v c="$c" '$1 == c && $3 == "exon" {print $1 "\t" $4 - 1 "\t" $5}' "$gtf" | sort -k2,2n \
    | awk 'BEGIN {OFS = "\t"} {if (NR > 1 && $2 <= e) {if ($3 > e) e = $3} else {if (NR > 1) print k, s, e; k = $1; s = $2; e = $3}} END {print k, s, e}' \
    > "$reg.tmp.$SLURM_ARRAY_TASK_ID"
  mv -n "$reg.tmp.$SLURM_ARRAY_TASK_ID" "$reg" 2>/dev/null || rm -f "$reg.tmp.$SLURM_ARRAY_TASK_ID"
fi
redi=$out/regions/$c.rediportal.tsv
if [ ! -s "$redi" ]; then
  zcat "$ext/rediportal/TABLE1_hg38_v3.txt.gz" | awk -F'\t' -v c="$c" 'NR > 1 && $2 == c {print $2 "\t" $3}' | sort -k2,2n -u > "$redi.tmp.$SLURM_ARRAY_TASK_ID"
  mv -n "$redi.tmp.$SLURM_ARRAY_TASK_ID" "$redi" 2>/dev/null || rm -f "$redi.tmp.$SLURM_ARRAY_TASK_ID"
fi
excl=$out/regions/$c.blacklist_mhc.tsv
if [ ! -s "$excl" ]; then
  { if [ -n "$bl" ]; then zcat -f "$bl" | awk -v c="$c" '$1 == c {print $1 "\t" $2 + 1 "\t" $3}'; fi
    if [ "$c" = chr6 ]; then printf 'chr6\t28500000\t33500000\n'; fi; } | sort -k2,2n > "$excl.tmp.$SLURM_ARRAY_TASK_ID"
  mv -n "$excl.tmp.$SLURM_ARRAY_TASK_ID" "$excl" 2>/dev/null || rm -f "$excl.tmp.$SLURM_ARRAY_TASK_ID"
fi

ls "$a1/bai/$cohort"/*.bam > "$out/calls/$cohort.$c.bams.txt"
bcftools mpileup -f "$ref" -R "$reg" -q 255 -Q 20 -d 1000 -a FORMAT/AD,FORMAT/DP -b "$out/calls/$cohort.$c.bams.txt" -Ou \
  | bcftools call -m -v -Ou \
  | bcftools view -v snps -m2 -M2 -i 'QUAL >= 30' -T "^$redi" -Ou \
  | { if [ -s "$excl" ]; then bcftools view -T "^$excl" -Ou; else cat; fi; } \
  | bcftools +setGT - -Ou -- -t q -n . -i 'FMT/DP < 10' \
  | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o "$out/calls/$cohort.$c.raw.vcf.gz"
# sample names are symlink paths "<dir>/<run>.bam" -> run IDs
bcftools query -l "$out/calls/$cohort.$c.raw.vcf.gz" | awk '{n = split($0, p, "/"); r = p[n]; sub(/\.bam$/, "", r); print $0 "\t" r}' \
  > "$out/calls/$cohort.$c.rename.tsv"
bcftools reheader -s "$out/calls/$cohort.$c.rename.tsv" -o "$out/calls/$cohort.$c.vcf.gz" "$out/calls/$cohort.$c.raw.vcf.gz"
rm "$out/calls/$cohort.$c.raw.vcf.gz"
bcftools index -t "$out/calls/$cohort.$c.vcf.gz"
echo "$cohort $c sites $(bcftools view -H "$out/calls/$cohort.$c.vcf.gz" | wc -l)"
