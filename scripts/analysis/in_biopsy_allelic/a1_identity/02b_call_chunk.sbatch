#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=24:00:00
#SBATCH --array=0-99
# A1 step 2 as an array: one cohort x one of 10 contiguous genome chunks of the
# identity targets per task, so BAM reads spread across nodes (one job on one node
# was seek-bound at ~75 MB/s). Same calling settings as 02_call_cohort.sh.
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%A_%a.out 02b_call_chunk.sbatch <targets_dir> <restricted_dir>
set -euo pipefail
tdir=$1; rs=$2
cohorts=(GSE126848 GSE130970 GSE135251 GSE162694 GSE167523 GSE174478 GSE193066 GSE213621 GSE240729 PRJNA512027)
cohort=${cohorts[$((SLURM_ARRAY_TASK_ID / 10))]}
k=$((SLURM_ARRAY_TASK_ID % 10))
ref=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa
module load bcftools/1.21 htslib/1.21 >/dev/null 2>&1
out=$rs/vcf_chunks; mkdir -p "$out"
n=$(wc -l < "$tdir/regions.tsv")
awk -v k="$k" -v n="$n" 'NR > int(k * n / 10) && NR <= int((k + 1) * n / 10)' "$tdir/regions.tsv" > "$out/$cohort.chunk$k.regions.tsv"
ls "$rs/bai/$cohort"/*.bam > "$out/$cohort.chunk$k.bams.txt"
bcftools mpileup -f "$ref" -R "$out/$cohort.chunk$k.regions.tsv" -q 255 -Q 20 -d 1000 \
    -a FORMAT/AD,FORMAT/DP -b "$out/$cohort.chunk$k.bams.txt" -Ou \
  | bcftools call -m -C alleles -T "$tdir/targets.tsv.gz" -Ou \
  | bcftools +setGT - -Oz -o "$out/$cohort.chunk$k.raw.vcf.gz" -- -t q -n . -i 'FMT/DP<10'
bcftools query -l "$out/$cohort.chunk$k.raw.vcf.gz" | awk '{n = split($0, p, "/"); r = p[n]; sub(/\.bam$/, "", r); print $0 "\t" r}' \
  > "$out/$cohort.chunk$k.rename.tsv"
bcftools reheader -s "$out/$cohort.chunk$k.rename.tsv" -o "$out/$cohort.chunk$k.vcf.gz" "$out/$cohort.chunk$k.raw.vcf.gz"
rm "$out/$cohort.chunk$k.raw.vcf.gz"
bcftools index -t "$out/$cohort.chunk$k.vcf.gz"
echo "$cohort chunk $k sites $(bcftools view -H "$out/$cohort.chunk$k.vcf.gz" | wc -l)"
