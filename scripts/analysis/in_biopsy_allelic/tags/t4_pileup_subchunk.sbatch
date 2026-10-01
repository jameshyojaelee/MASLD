#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=48:00:00
#SBATCH --array=0-19
# T4 rescue: one slow cohort chunk split into 20 contiguous sub-chunks, same settings as
# t4_pileup.sbatch (MAPQ 255, BQ 20, 10,000 reads per library per site, BAQ default).
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%A_%a.out t4_pileup_subchunk.sbatch <ad_dir> <cohort> <chunk>
set -euo pipefail
ad=$1; cohort=$2; k=$3; s=$SLURM_ARRAY_TASK_ID
ref=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa
module load bcftools/1.21 htslib/1.21 >/dev/null 2>&1
reg=$ad/$cohort.chunk$k.regions.tsv
n=$(wc -l < "$reg")
mkdir -p "$ad/sub"
awk -v s="$s" -v n="$n" 'NR > int(s * n / 20) && NR <= int((s + 1) * n / 20)' "$reg" > "$ad/sub/$cohort.chunk$k.sub$s.regions.tsv"
[ -s "$ad/sub/$cohort.chunk$k.sub$s.regions.tsv" ] || { echo "empty sub-chunk"; exit 0; }
bcftools mpileup -f "$ref" -R "$ad/sub/$cohort.chunk$k.sub$s.regions.tsv" -q 255 -Q 20 -d 10000 -a FORMAT/AD \
    -b "$ad/$cohort.chunk$k.bams.txt" -Ou \
  | bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' | gzip > "$ad/sub/$cohort.chunk$k.sub$s.ad.tsv.gz"
echo "$cohort chunk $k sub $s rows $(zcat "$ad/sub/$cohort.chunk$k.sub$s.ad.tsv.gz" | wc -l) sites $(wc -l < "$ad/sub/$cohort.chunk$k.sub$s.regions.tsv")"
