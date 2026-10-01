#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
# A9: allele depths at six named coding / risk SNVs (coding_variants.tsv; positions from
# dbSNP b157 within each gene's GENCODE v49 span, alleles from the 1kGP panel) on every
# library, same settings as T4 (MAPQ 255, BQ 20, 10,000 reads per library per site).
# Usage: sbatch --job-name=<n> --output=<log> a9_pileup.sbatch <a1_restricted_dir> <out_restricted_dir>
set -euo pipefail
a1=$1; out=$2
here=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/in_biopsy_allelic/a9_coding
ref=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa
module load bcftools/1.21 htslib/1.21 >/dev/null 2>&1
mkdir -p "$out"
awk -F'\t' 'NR > 1 {print $3 "\t" $4 "\t" $4}' "$here/coding_variants.tsv" | sort -k1,1V -k2,2n > "$out/regions.tsv"
for c in GSE126848 GSE130970 GSE135251 GSE162694 GSE167523 GSE174478 GSE193066 GSE213621 GSE240729 PRJNA512027; do
  ls "$a1/bai/$c"/*.bam > "$out/$c.bams.txt"
  bcftools mpileup -f "$ref" -R "$out/regions.tsv" -q 255 -Q 20 -d 10000 -a FORMAT/AD -b "$out/$c.bams.txt" -Ou \
    | bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' | gzip > "$out/$c.ad.tsv.gz"
  echo "$c rows $(zcat "$out/$c.ad.tsv.gz" | wc -l)"
done
