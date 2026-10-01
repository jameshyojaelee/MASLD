#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --array=0-47
# A1 step 2, sub-chunks for slow cohort chunks: each listed chunk's regions file is
# split into 8 contiguous parts; same calling settings as 02b_call_chunk.sbatch.
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%A_%a.out 02c_call_subchunk.sbatch <targets_dir> <restricted_dir> <cohort> "<chunk list>"
set -euo pipefail
tdir=$1; rs=$2; cohort=$3; read -r -a chunks <<< "$4"
k=${chunks[$((SLURM_ARRAY_TASK_ID / 8))]}
s=$((SLURM_ARRAY_TASK_ID % 8))
ref=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa
module load bcftools/1.21 htslib/1.21 >/dev/null 2>&1
out=$rs/vcf_chunks; sub=$out/sub; mkdir -p "$sub"
reg=$out/$cohort.chunk$k.regions.tsv
n=$(wc -l < "$reg")
awk -v s="$s" -v n="$n" 'NR > int(s * n / 8) && NR <= int((s + 1) * n / 8)' "$reg" > "$sub/$cohort.chunk$k.sub$s.regions.tsv"
bcftools mpileup -f "$ref" -R "$sub/$cohort.chunk$k.sub$s.regions.tsv" -q 255 -Q 20 -d 1000 \
    -a FORMAT/AD,FORMAT/DP -b "$out/$cohort.chunk$k.bams.txt" -Ou \
  | bcftools call -m -C alleles -T "$tdir/targets.tsv.gz" -Ou \
  | bcftools +setGT - -Oz -o "$sub/$cohort.chunk$k.sub$s.vcf.gz" -- -t q -n . -i 'FMT/DP<10'
echo "$cohort chunk $k sub $s sites $(bcftools view -H "$sub/$cohort.chunk$k.sub$s.vcf.gz" | wc -l)"
