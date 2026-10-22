#!/bin/bash
# CGB4: MAGMA gene-property (M1) and conditional system gene-set (M2) analyses.
# Args: CG0 out dir (genes.tsv, gene_system_weights.tsv.gz), CG0 MAGMA run dir, output dir.
set -euo pipefail
CG0=$1; MRUN=$2; OUT=$3
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
HERE=$(dirname "$(readlink -f "$0")")
MAGMA=$ROOT/tools/magma/magma
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1
mkdir -p $OUT
$PY $HERE/cgb4_build_covar.py --cg0 $CG0 --out $OUT
for t in ukbb_alt ukbb_ast mvp_alt_eur mvp_ast_eur pazoki_pdff ghouse_cirrhosis yengo_bmi; do
  R=$MRUN/$t.genes.raw
  test -s $R || { echo "missing $R"; continue; }
  $MAGMA --gene-results $R --gene-covar $OUT/covar.txt \
    --model direction=pos condition-hide=P_dis,AveExpr,LOEUF analyse=list,P_cont --out $OUT/${t}.M1_cont > /dev/null &
  $MAGMA --gene-results $R --gene-covar $OUT/covar.txt \
    --model direction=pos condition-hide=P_cont,AveExpr,LOEUF analyse=list,P_dis --out $OUT/${t}.M1_dis > /dev/null &
  $MAGMA --gene-results $R --gene-covar $OUT/covar.txt --set-annot $OUT/systems.sets \
    --model direction=pos condition-hide=P_cont,P_dis,AveExpr,LOEUF analyse=sets --out $OUT/${t}.M2_sys > /dev/null &
  wait
done
$PY $HERE/cgb4_collate.py --out $OUT
