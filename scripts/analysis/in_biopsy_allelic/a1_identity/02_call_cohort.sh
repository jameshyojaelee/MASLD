#!/usr/bin/env bash
# A1 step 2: call the identity-panel SNVs in one cohort's existing STAR BAMs.
# Usage: 02_call_cohort.sh <cohort> <targets_dir> <out_dir> <index_dir> [max_bams]
# The existing BAMs are mostly unindexed. Indexes are written to <index_dir>
# beside symlinks to the BAMs, never next to the original BAMs.
# Unique STAR alignments only (MAPQ 255); base quality >= 20; depth cap 1,000
# per file; genotypes with DP < 10 set to missing. Reads no outcome.
set -euo pipefail

cohort=$1
tdir=$2
out=$3
idx=$4
max_bams=${5:-0}
root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/results
ref=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa
samtools=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/samtools

mkdir -p "$out" "$idx/$cohort"
ls "$root/$cohort"/alignments/star/*/*.Aligned.sortedByCoord.out.bam > "$out/$cohort.bampaths.txt"
if [ "$max_bams" -gt 0 ]; then head -n "$max_bams" "$out/$cohort.bampaths.txt" > "$out/$cohort.tmp" && mv "$out/$cohort.tmp" "$out/$cohort.bampaths.txt"; fi
test -s "$out/$cohort.bampaths.txt"

# bcftools rejects "##idx##" in a BAM list, so each BAM is reached through a
# symlink in the index dir with its index beside it (originals untouched).
while read -r bam; do
  run=$(basename "$bam" .Aligned.sortedByCoord.out.bam)
  link="$idx/$cohort/$run.bam"
  ln -sfn "$bam" "$link"
  [ -s "$link.bai" ] || "$samtools" index -o "$link.bai" "$link"
  printf '%s\t%s\n' "$link" "$run"
done < "$out/$cohort.bampaths.txt" > "$out/$cohort.bams_idx_run.tsv"
cut -f1 "$out/$cohort.bams_idx_run.tsv" > "$out/$cohort.bams.txt"

bcftools mpileup -f "$ref" -R "$tdir/regions.tsv" -q 255 -Q 20 -d 1000 \
    -a FORMAT/AD,FORMAT/DP -b "$out/$cohort.bams.txt" -Ou \
  | bcftools call -m -C alleles -T "$tdir/targets.tsv.gz" -Ou \
  | bcftools +setGT - -Oz -o "$out/$cohort.raw.vcf.gz" -- -t q -n . -i 'FMT/DP<10'

# Sample names are the BAM paths as bcftools saw them; map them to run IDs.
bcftools query -l "$out/$cohort.raw.vcf.gz" > "$out/$cohort.rawnames.txt"
awk -F'\t' 'NR==FNR {m[$1]=$2; next}
            {print $0"\t"(($0 in m) ? m[$0] : "UNMAPPED")}' \
    "$out/$cohort.bams_idx_run.tsv" "$out/$cohort.rawnames.txt" > "$out/$cohort.rename.tsv"
if grep -q UNMAPPED "$out/$cohort.rename.tsv"; then echo "unmapped sample names in $cohort" >&2; exit 1; fi
bcftools reheader -s "$out/$cohort.rename.tsv" -o "$out/$cohort.vcf.gz" "$out/$cohort.raw.vcf.gz"
rm "$out/$cohort.raw.vcf.gz"
bcftools index -t "$out/$cohort.vcf.gz"
echo "$cohort samples=$(bcftools query -l "$out/$cohort.vcf.gz" | wc -l) sites=$(bcftools view -H "$out/$cohort.vcf.gz" | wc -l)"
