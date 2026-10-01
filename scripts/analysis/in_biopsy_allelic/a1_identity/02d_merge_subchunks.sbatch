#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
# Concatenate the 8 sub-chunks of each listed chunk into the chunk VCF that
# 03b_identity.sbatch expects, with run IDs as sample names (as in 02b).
# Usage: sbatch --job-name=<n> --dependency=afterok:<02c> --output=<log> 02d_merge_subchunks.sbatch <restricted_dir> <cohort> "<chunk list>"
set -euo pipefail
rs=$1; cohort=$2; read -r -a chunks <<< "$3"
module load bcftools/1.21 htslib/1.21 >/dev/null 2>&1
out=$rs/vcf_chunks; sub=$out/sub
for k in "${chunks[@]}"; do
  bcftools concat -Oz -o "$sub/$cohort.chunk$k.raw.vcf.gz" $(for s in 0 1 2 3 4 5 6 7; do printf '%s ' "$sub/$cohort.chunk$k.sub$s.vcf.gz"; done)
  bcftools query -l "$sub/$cohort.chunk$k.raw.vcf.gz" | awk '{n = split($0, p, "/"); r = p[n]; sub(/\.bam$/, "", r); print $0 "\t" r}' \
    > "$out/$cohort.chunk$k.rename.tsv"
  bcftools reheader -s "$out/$cohort.chunk$k.rename.tsv" -o "$out/$cohort.chunk$k.vcf.gz" "$sub/$cohort.chunk$k.raw.vcf.gz"
  bcftools index -t "$out/$cohort.chunk$k.vcf.gz"
  echo "$cohort chunk $k sites $(bcftools view -H "$out/$cohort.chunk$k.vcf.gz" | wc -l)"
done
