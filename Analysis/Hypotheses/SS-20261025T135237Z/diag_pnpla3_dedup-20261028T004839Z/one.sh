#!/bin/bash
# Post-hoc diagnostic (not prespecified): allele counts at 4 sites with all reads vs after duplicate removal.
set -euo pipefail
S=$1; O=$2; ST=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/samtools
B=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/Human/Patient_Cohorts/results/SteatoSITE/alignments/star/$S/$S.Aligned.sortedByCoord.out.bam
FA=/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49/genome.primary_chr.fa
t=$(mktemp -d -p $O/tmp); trap 'rm -rf $t' EXIT
$ST view -b -q 20 -F 0x904 $B chr22:43918847-43938847 chr2:27498073-27518073 chr19:44898684-44918684 chr19:19258740-19278740 > $t/r.bam
$ST collate -O -u $t/r.bam $t/c | $ST fixmate -m -u - - | $ST sort -u -o $t/f.bam -
$ST markdup -r -s $t/f.bam $t/d.bam 2> $t/md.txt; $ST index $t/f.bam; $ST index $t/d.bam
for k in all:f dedup:d; do
  $ST mpileup -f $FA -l $O/sites.tsv -q 20 -Q 20 -B -d 100000 $t/${k#*:}.bam 2>/dev/null | awk -v s=$S -v m=${k%:*} 'BEGIN{OFS="\t"}{b=toupper($5); gsub(/\^./,"",b); gsub(/[$]/,"",b); r=gsub(/[.,]/,"",b); a=length(b)-gsub(/[*]/,"",b); print s,m,$1,$2,$4,r,b}'
done
