#!/bin/bash
#SBATCH --job-name=DIG1magma69
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/logs/%j.out
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RUN=$ROOT/Analysis/Hypotheses/DIG1-20261024T113951Z
OUT=$RUN/out_state_magma
MRUN=$ROOT/Analysis/Hypotheses/CG0-20261022T122300Z/magma/run-20261023T022444Z
COV=$ROOT/Analysis/Hypotheses/CGB4-20261022T122538Z/out-20261023T032051Z/covar.txt
MAGMA=$ROOT/tools/magma/magma
sha256sum $RUN/scripts/69_state_magma_sets.py $COV > $RUN/script_sha256_job69.txt
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -I $RUN/scripts/69_state_magma_sets.py $OUT
for t in ukbb_alt ukbb_ast mvp_alt_eur mvp_ast_eur pazoki_pdff ghouse_cirrhosis yengo_bmi; do
  $MAGMA --gene-results $MRUN/$t.genes.raw --gene-covar $COV --set-annot $OUT/state.sets \
    --model direction=pos condition-hide=P_cont,P_dis,AveExpr,LOEUF analyse=sets --out $OUT/$t.state > /dev/null &
  $MAGMA --gene-results $MRUN/$t.genes.raw --gene-covar $COV --set-annot $OUT/state.sets \
    --model direction=pos condition-hide=AveExpr,LOEUF analyse=sets --out $OUT/$t.state_noPdis > /dev/null &
  # sensitivity: drop genes with gene P < 2.5e-6 in this trait
  awk 'NR>1 && $9 < 2.5e-6 {print $1}' $MRUN/$t.genes.out > $OUT/$t.gwsig_genes.txt
  awk 'NR==FNR{d[$1]=1; next} {o=$1; for(i=2;i<=NF;i++) if(!($i in d)) o=o" "$i; print o}' $OUT/$t.gwsig_genes.txt $OUT/state.sets > $OUT/$t.state_nosig.sets
  $MAGMA --gene-results $MRUN/$t.genes.raw --gene-covar $COV --set-annot $OUT/$t.state_nosig.sets \
    --model direction=pos condition-hide=P_cont,P_dis,AveExpr,LOEUF analyse=sets --out $OUT/$t.state_nosig > /dev/null &
  wait
done
grep -H -v "^#" $OUT/*.gsa.out | grep -v VARIABLE
