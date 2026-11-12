#!/bin/bash
# GSET8 follow-ups, post hoc except where noted in RESULT_GSET8.md.
# (a) locus exclusion and top-gene removal; (b) MVP cirrhosis (phecode 571.51) replication;
# (d) hepatocyte-tested background. Args: primary out dir, new follow-up out dir.
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
HERE=$(dirname "$(readlink -f "$0")")
PRIM=$1; OUT=$2
cd $ROOT
test ! -e $OUT || { echo "output exists: $OUT"; exit 1; }
mkdir -p $OUT
MAGMA=$ROOT/tools/magma/magma
REF=$ROOT/tools/magma/ref/g1000_eur
ANNOT=$ROOT/tools/magma/work/pops/magma_annot_ensgid.genes.annot
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1
MVP=GWAS/MR_Data/MVP/raw/phecodes/MVP_R4.1000G_AGR.Phe_571_51.EUR.GIA.dbGaP.txt.gz
STOUF=Analysis/Hypotheses/CONF1-20261109T221846Z/out_c3b_ambient/hep_corrected_stouffer.tsv.gz
{ $MAGMA --version 2>&1 | head -1; hostname; date -u; } > $OUT/versions.txt
sha256sum $MVP $STOUF $PRIM/*.genes.raw $PRIM/sets_primary.txt $PRIM/sets_with_random.txt > $OUT/input_sha256.txt
sha256sum $HERE/*.sh $HERE/*.py > $OUT/script_sha256.txt

# (b) MVP input: rs IDs, imputation r2 >= 0.9 (ENZFAT1 rule), 0 < P <= 1 (CG0 rule), first row per rsID, N = num_samples
zcat $MVP | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; print "SNP\tP\tN"; next}
  {tot++; s=$h["SNP_ID"]; p=$h["pval"]; r=$h["r2"];
   if(s !~ /^rs[0-9]+$/){bid++; next} if(!(r+0>=0.9)){br2++; next}
   if(!(p+0>0 && p+0<=1)){bp++; next} if(s in seen){dup++; next}
   seen[s]=1; k++; print s"\t"p"\t"$h["num_samples"]}
  END{printf "mvp_cirrhosis rows=%d kept=%d dropped_id=%d dropped_r2=%d dropped_p=%d dropped_dup=%d\n",tot,k,bid,br2,bp,dup > "/dev/stderr"}' \
  > $OUT/mvp_cirrhosis.magma_input.tsv 2> $OUT/input_filter_counts.txt
cat $OUT/input_filter_counts.txt
$MAGMA --bfile $REF --pval $OUT/mvp_cirrhosis.magma_input.tsv use=SNP,P ncol=N \
  --gene-annot $ANNOT --out $OUT/mvp_cirrhosis > $OUT/mvp_cirrhosis.gene_stdout 2>&1

$PY $HERE/gset8_followup_lists.py --prim $PRIM --out $OUT

raw() { if [ $1 = mvp_cirrhosis ]; then echo $OUT/mvp_cirrhosis.genes.raw; else echo $PRIM/$1.genes.raw; fi; }
run() {  # tag trait extra-args...
  local tag=$1 t=$2; shift 2
  $MAGMA --gene-results $(raw $t) "$@" --out $OUT/$t.$tag > $OUT/$t.$tag.stdout 2>&1
}
for t in pdff cirrhosis mvp_cirrhosis; do
  run excl_loci $t --set-annot $PRIM/sets_primary.txt --settings gene-exclude=$OUT/exclude_loci_genes.txt --model direction=pos &
  run heptested $t --set-annot $PRIM/sets_primary.txt --settings gene-include=$OUT/include_heptested.txt --model direction=pos &
  run heptested_aveexpr $t --gene-covar $PRIM/covar_aveexpr.txt --set-annot $PRIM/sets_primary.txt \
    --settings gene-include=$OUT/include_heptested.txt --model direction=pos condition-hide=AveExpr analyse=sets &
  run heptested_down $t --set-annot $PRIM/sets_primary.txt --settings gene-include=$OUT/include_heptested_down.txt --model direction=pos &
  wait
done
run drop_top cirrhosis --set-annot $OUT/sets_drop_top.txt --model direction=pos &
run drop_top mvp_cirrhosis --set-annot $OUT/sets_drop_top.txt --model direction=pos &
run primary mvp_cirrhosis --set-annot $PRIM/sets_with_random.txt --model direction=pos &
run sens_aveexpr mvp_cirrhosis --gene-covar $PRIM/covar_aveexpr.txt --set-annot $PRIM/sets_primary.txt \
  --model direction=pos condition-hide=AveExpr analyse=sets &
wait
$PY $HERE/gset8_followup_collate.py --out $OUT
