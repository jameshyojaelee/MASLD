#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=48:00:00
#SBATCH --array=0-199
# T4: allele depths at the final tag SNVs, one cohort x one of 20 genome chunks per
# task, on the existing STAR BAMs through the A1 symlinked indexes. MAPQ 255, BQ 20,
# depth cap 10,000 reads per library per site (position order, allele-independent).
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%A_%a.out t4_pileup.sbatch <tags_final.tsv> <a1_restricted_dir> <out_restricted_dir>
set -euo pipefail
tags=$1; a1=$2; out=$3
cohorts=(GSE126848 GSE130970 GSE135251 GSE162694 GSE167523 GSE174478 GSE193066 GSE213621 GSE240729 PRJNA512027)
cohort=${cohorts[$((SLURM_ARRAY_TASK_ID / 20))]}
k=$((SLURM_ARRAY_TASK_ID % 20))
ref=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa
module load bcftools/1.21 htslib/1.21 >/dev/null 2>&1
mkdir -p "$out/ad"
sites=$out/ad/sites.sorted.tsv
if [ ! -s "$sites" ]; then
  awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} {print $h["chrom"] "\t" $h["pos"] "\t" $h["pos"]}' "$tags" \
    | sort -k1,1V -k2,2n -u > "$sites.tmp.$SLURM_ARRAY_TASK_ID"
  mv -n "$sites.tmp.$SLURM_ARRAY_TASK_ID" "$sites" 2>/dev/null || rm -f "$sites.tmp.$SLURM_ARRAY_TASK_ID"
fi
n=$(wc -l < "$sites")
awk -v k="$k" -v n="$n" 'NR > int(k * n / 20) && NR <= int((k + 1) * n / 20)' "$sites" > "$out/ad/$cohort.chunk$k.regions.tsv"
[ -s "$out/ad/$cohort.chunk$k.regions.tsv" ] || { echo "empty chunk"; exit 0; }
ls "$a1/bai/$cohort"/*.bam > "$out/ad/$cohort.chunk$k.bams.txt"
bcftools mpileup -f "$ref" -R "$out/ad/$cohort.chunk$k.regions.tsv" -q 255 -Q 20 -d 10000 -a FORMAT/AD \
    -b "$out/ad/$cohort.chunk$k.bams.txt" -Ou \
  | bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' | gzip > "$out/ad/$cohort.chunk$k.ad.tsv.gz"
echo "$cohort chunk $k rows $(zcat "$out/ad/$cohort.chunk$k.ad.tsv.gz" | wc -l)"
