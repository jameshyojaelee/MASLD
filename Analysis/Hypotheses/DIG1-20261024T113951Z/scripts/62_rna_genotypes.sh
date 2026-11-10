#!/bin/bash
# DIG1 job 62. Genotypes at 10 liver-trait variants called from the bulk RNA-seq reads of the 774 release
# samples (STAR 2.7.11b BAMs, GENCODE v49 / GRCh38). Positions and alleles from Ensembl REST/VEP
# (downloads_untrusted/ensembl_variant_positions, SHA256SUMS) and checked against the reference FASTA.
# Step 1: a BAI index per BAM, written to <out>/bai (source BAMs are only read).
# Step 2: per sample, bcftools mpileup at the targets (MAPQ >= 20, i.e. STAR unique; base Q >= 20; no BAQ;
#         max depth 100000) and bcftools call -m; FORMAT/AD kept for count-based calls (job 62b).
# Usage: 62_rna_genotypes.sh <new_out_dir>
set -euo pipefail
OUT=$1; [ -e "$OUT" ] && { echo "refusing to overwrite $OUT"; exit 1; }
mkdir -p "$OUT/bai" "$OUT/vcf" "$OUT/logs"
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z
B=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/Human/Patient_Cohorts/results
REL=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls
FA=/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49/genome.primary_chr.fa
S=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/samtools
BCF=/gpfs/commons/home/jameslee/micromamba/envs/allelic_qtl/bin/bcftools
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript -e "m <- readRDS('$REL/meta_matched.rds'); write.table(m[, c('sample_id', 'dataset')], '$OUT/samples.tsv', sep = '\t', quote = FALSE, row.names = FALSE)"
cp "$RUN/scripts/62_targets.tsv" "$RUN/scripts/62_targets_annotated.tsv" "$OUT/"
one() {
  s=$1; d=$2
  bam=$B/$d/alignments/star/$s/$s.Aligned.sortedByCoord.out.bam
  [ -s "$OUT/bai/$s.bai" ] || $S index -@ 1 -o "$OUT/bai/$s.bai" "$bam"
  $BCF mpileup -f "$FA" -R "$OUT/62_targets.tsv" -a FORMAT/AD,FORMAT/DP -q 20 -Q 20 -B -d 100000 "$bam##idx##$OUT/bai/$s.bai" 2> "$OUT/logs/$s.mpileup.err" \
    | $BCF call -m -Oz -o "$OUT/vcf/$s.vcf.gz" 2> "$OUT/logs/$s.call.err"
  echo "$s done"
}
export -f one; export OUT B FA S BCF
tail -n +2 "$OUT/samples.tsv" | xargs -P 8 -L 1 bash -c 'one $0 $1'
n=$(ls "$OUT"/vcf/*.vcf.gz | wc -l); echo "vcf files: $n"
