#!/bin/bash
#SBATCH --job-name=DIG1null72
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/logs/%j.out
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z
OUT=$RUN/out_ase_null_sites; [ -e $OUT ] && exit 1; mkdir -p $OUT/vcf $OUT/logs
PY="/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -I"
BCF=/gpfs/commons/home/jameslee/micromamba/envs/allelic_qtl/bin/bcftools
K=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external/allelic_refs/kgp_highcov
sha256sum $RUN/scripts/72a_null_utr_sites.py $RUN/scripts/72b_null_ase.R > $RUN/script_sha256_job72.txt
$PY $RUN/scripts/72a_null_utr_sites.py regions $OUT/regions.tsv
cut -f1-3 $OUT/regions.tsv > $OUT/regions.bed3
for c in $(seq 1 22); do
  $BCF view -v snps -m2 -M2 -R $OUT/regions.bed3 $K/1kGP_high_coverage_Illumina.chr$c.filtered.SNV_INDEL_SV_phased_panel.vcf.gz 2>/dev/null \
    | $BCF query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF_EUR\n'
done > $OUT/kgp_snvs.tsv
$PY $RUN/scripts/72a_null_utr_sites.py pick $OUT/kgp_snvs.tsv $OUT/regions.tsv $OUT/targets_annot.tsv
cut -f1-3 $OUT/targets_annot.tsv | sort -k1,1V -k2,2n > $OUT/targets.tsv
B=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/Human/Patient_Cohorts/results
FA=/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49/genome.primary_chr.fa
one() {
  s=$1; d=$2; bam=$B/$d/alignments/star/$s/$s.Aligned.sortedByCoord.out.bam
  $BCF mpileup -f "$FA" -R "$OUT/targets.tsv" -a FORMAT/AD,FORMAT/DP -q 20 -Q 20 -B -d 100000 "$bam##idx##$RUN/out_rna_genotypes/bai/$s.bai" 2> "$OUT/logs/$s.err" \
    | $BCF query -f "$s\t%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n" > "$OUT/vcf/$s.tsv"
}
export -f one; export OUT B FA BCF RUN
tail -n +2 $RUN/out_rna_genotypes/samples.tsv | xargs -P 8 -L 1 bash -c 'one $0 $1'
cat $OUT/vcf/*.tsv > $OUT/allele_counts.tsv; wc -l $OUT/allele_counts.tsv
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $RUN/scripts/72b_null_ase.R $OUT $RUN/out_panx1_ase/utr_het_counts.tsv
