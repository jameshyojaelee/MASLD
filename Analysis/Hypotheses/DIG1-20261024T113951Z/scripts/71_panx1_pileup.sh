#!/bin/bash
# DIG1 job 71a. Allele counts at 7 PANX1 SNPs (GRCh38) in the 774 release bulk RNA BAMs, for allele-specific
# expression. SNPs: exon 1 rs? 94128854 C>T (GTEx liver eQTL slope +0.21, separate weak signal; r -0.26 with the
# coloc haplotype); intronic coloc-haplotype SNPs 94132981 T>G, 94137172 C>T, 94144210 A>G (Broadaway eQTL lead,
# rs4593973), 94150836 C>G (GTEx lead, rs12273690); 3'UTR tags 94181312 G>C and 94181824 T>G (EUR r +0.74 with the
# lead, 1kGP high-coverage phased panel, 633 EUR). Positions hg19 -> hg38 by a +266,834 offset checked on rs12273690,
# then confirmed in the 1kGP hg38 VCF. Same mpileup settings as job 62 (MAPQ >= 20, BQ >= 20, no BAQ), BAI from job 62.
# Usage: 71_panx1_pileup.sh <new_out_dir>
set -euo pipefail
OUT=$1; [ -e "$OUT" ] && { echo "refusing to overwrite $OUT"; exit 1; }
mkdir -p "$OUT/vcf" "$OUT/logs"
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z
B=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/Human/Patient_Cohorts/results
FA=/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49/genome.primary_chr.fa
BCF=/gpfs/commons/home/jameslee/micromamba/envs/allelic_qtl/bin/bcftools
SAM=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/samtools
cp "$RUN/scripts/71_targets.tsv" "$OUT/"
while read c p e; do echo -n "$c:$p ref "; $SAM faidx "$FA" "$c:$p-$p" | tail -1; done < "$OUT/71_targets.tsv" > "$OUT/ref_check.txt"
one() {
  s=$1; d=$2
  bam=$B/$d/alignments/star/$s/$s.Aligned.sortedByCoord.out.bam
  $BCF mpileup -f "$FA" -R "$OUT/71_targets.tsv" -a FORMAT/AD,FORMAT/DP -q 20 -Q 20 -B -d 100000 "$bam##idx##$RUN/out_rna_genotypes/bai/$s.bai" 2> "$OUT/logs/$s.err" \
    | $BCF query -f "$s\t%POS\t%REF\t%ALT\t[%AD]\n" > "$OUT/vcf/$s.tsv"
}
export -f one; export OUT B FA BCF RUN
tail -n +2 "$RUN/out_rna_genotypes/samples.tsv" | xargs -P 8 -L 1 bash -c 'one $0 $1'
cat "$OUT"/vcf/*.tsv > "$OUT/allele_counts.tsv"; wc -l "$OUT/allele_counts.tsv"; cat "$OUT/ref_check.txt"
