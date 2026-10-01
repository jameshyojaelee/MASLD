#!/bin/bash
#SBATCH --partition=io
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=72:00:00
# Same downloads as download_refs.sbatch, 8 wget streams in parallel; wget -c
# resumes files partly fetched by the serial job (0903, cancelled for speed).
# Usage: sbatch --job-name=<n> --output=<log>/%j.out download_refs_parallel.sbatch
set -euo pipefail
ext=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external
refs=$ext/allelic_refs
kgp=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV
mkdir -p "$refs/kgp_highcov" "$refs/rediportal" "$refs/gtex_v8" "$ext/geuvadis/fastq"

jobs=$ext/download_jobs_$(date -u +%Y%m%dT%H%M%SZ).tsv
{
  for c in $(seq 1 22); do
    f=1kGP_high_coverage_Illumina.chr$c.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
    printf '%s\t%s\n' "$refs/kgp_highcov" "$kgp/$f" "$refs/kgp_highcov" "$kgp/$f.tbi"
  done
  f=1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz
  printf '%s\t%s\n' "$refs/kgp_highcov" "$kgp/$f" "$refs/kgp_highcov" "$kgp/$f.tbi"
  printf '%s\t%s\n' "$refs/rediportal" http://rediportal.cloud.ba.infn.it/download/TABLE1_hg38_v3.txt.gz
  printf '%s\t%s\n' "$refs/gtex_v8" https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar
  awk -F'\t' 'NR > 1 {print $5}' "$ext/geuvadis/selected_30_eur.tsv" | tr ';' '\n' \
    | awk -v d="$ext/geuvadis/fastq" '{print d "\thttp://" $0}'
} > "$jobs"
echo "$(wc -l < "$jobs") files"

awk -F'\t' '{print $1 "\t" $2}' "$jobs" | xargs -P 8 -L 1 sh -c 'cd "$0" && wget -c -q "$1" || echo "FAILED $1"'

cd "$refs/gtex_v8"
tar -tf GTEx_Analysis_v8_eQTL.tar | grep -i "Liver" > liver_members.txt
tar -xf GTEx_Analysis_v8_eQTL.tar -T liver_members.txt
for d in "$refs/kgp_highcov" "$refs/rediportal" "$refs/gtex_v8" "$ext/geuvadis/fastq"; do
  (cd "$d" && find . -type f ! -name 'SHA256SUMS*' ! -name '*.txt' -print0 | xargs -0 sha256sum > "SHA256SUMS_$(date -u +%Y%m%dT%H%M%SZ).txt")
done
echo downloads done
