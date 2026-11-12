#!/bin/bash
# GSET8: MAGMA gene analysis (CG0 settings) for PDFF and Ghouse cirrhosis, then competitive
# gene-set tests of hep_up, hep_down, complexI with size-matched random sets.
# Arg: output directory (new).
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
HERE=$(dirname "$(readlink -f "$0")")
OUT=$1
cd $ROOT
test ! -e $OUT || { echo "output exists: $OUT"; exit 1; }
mkdir -p $OUT
MAGMA=$ROOT/tools/magma/magma
REF=$ROOT/tools/magma/ref/g1000_eur
ANNOT=$ROOT/tools/magma/work/pops/magma_annot_ensgid.genes.annot
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1
PDFF=GWAS/MR_Data/GCST90267352_PDFF_Pazoki2022.tsv.gz
CIRR=GWAS/MR_Data/Ghouse_Cirrhosis/GCST90319877_harmonised_hg38.tsv.gz
PART=Analysis/Hypotheses/CONF1-20261109T221846Z/out_c3b_ambient/final_partition_posthoc.tsv.gz
SYM=tools/pops/example/data/utils/gene_annot_jun10.txt
CG0G=Analysis/Hypotheses/CG0-20261022T122300Z/out-20261022T122649Z/genes.tsv

{ $MAGMA --version 2>&1 | head -1; awk --version | head -1; $PY --version; hostname; date -u; } > $OUT/versions.txt
$PY -m pip freeze > $OUT/pip_freeze.txt
sha256sum $PDFF $CIRR $PART $SYM $CG0G $ANNOT $REF.bed $REF.bim $REF.fam $MAGMA > $OUT/input_sha256.txt
sha256sum $HERE/*.sh $HERE/*.py $(dirname $HERE)/PRESPEC_GSET8.md > $OUT/script_sha256.txt

# SNP P N input. Header-based columns; 0 < P <= 1 (CG0 rule); rs IDs only; first row per rsID.
make_input() {  # file snpcol pcol ncol out
  zcat $1 | awk -F'\t' -v S=$2 -v P=$3 -v N=$4 -v TAG=$(basename $5) '
    NR==1{for(i=1;i<=NF;i++)h[$i]=i; print "SNP\tP\tN"; next}
    {tot++; s=$h[S]; p=$h[P]; n=$h[N];
     if(s !~ /^rs[0-9]+$/){bad_id++; next}
     if(!(p+0>0 && p+0<=1)){bad_p++; next}
     if(s in seen){dup++; next}
     seen[s]=1; kept++; print s"\t"p"\t"n}
    END{printf "%s rows=%d kept=%d dropped_id=%d dropped_p=%d dropped_dup=%d\n", TAG, tot, kept, bad_id, bad_p, dup > "/dev/stderr"}' > $5
}
make_input $PDFF variant_id p_value n $OUT/pdff.magma_input.tsv 2>> $OUT/input_filter_counts.txt &
make_input $CIRR rsid p_value samplesize $OUT/cirrhosis.magma_input.tsv 2>> $OUT/input_filter_counts.txt &
wait
cat $OUT/input_filter_counts.txt

for t in pdff cirrhosis; do
  $MAGMA --bfile $REF --pval $OUT/$t.magma_input.tsv use=SNP,P ncol=N \
    --gene-annot $ANNOT --out $OUT/$t > $OUT/$t.gene_stdout 2>&1 &
done
wait
grep -H -E "valid SNP|genes" $OUT/pdff.log $OUT/cirrhosis.log | head -20

$PY $HERE/gset8_build_sets.py --out $OUT

for t in pdff cirrhosis; do
  $MAGMA --gene-results $OUT/$t.genes.raw --set-annot $OUT/sets_with_random.txt \
    --model direction=pos --out $OUT/$t.primary > $OUT/$t.primary_stdout 2>&1 &
  $MAGMA --gene-results $OUT/$t.genes.raw --gene-covar $OUT/covar_aveexpr.txt --set-annot $OUT/sets_primary.txt \
    --model direction=pos condition-hide=AveExpr analyse=sets --out $OUT/$t.sens_aveexpr > $OUT/$t.sens_stdout 2>&1 &
done
wait

$PY $HERE/gset8_collate.py --out $OUT
