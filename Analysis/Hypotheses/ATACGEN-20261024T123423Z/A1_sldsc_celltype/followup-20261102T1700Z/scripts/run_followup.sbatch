#!/bin/bash
#SBATCH --job-name=A1_fu
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=24:00:00
# A1 follow-up (PRESPEC.md): Roadmap annotations (01), LD scores, model files (02),
# paired AST/ALT sumstats (03), regressions (04), summary (05).
# Usage: sbatch --output=<F>/logs/%x_%j.out --error=<F>/logs/%x_%j.err run_followup.sbatch <F>
set -euo pipefail
F="${1:?follow-up dir required}"
test -f "$F/PRESPEC.md"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
A1=$(dirname "$F")
H22=$ROOT/Analysis/Hypotheses/H22-20261009T015913Z
RENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
LENV=/gpfs/commons/home/jameslee/micromamba/envs/ldsc
REF=$ROOT/GWAS/ldsc_ref
PLINK=$ROOT/data/gsmap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC.
FRQ=$REF/sldsc/1000G_Phase3_frq/1000G.EUR.QC.
BASE=$REF/sldsc/baselineLD_v2.2/baselineLD.
export F A1 H22 RENV LENV PLINK FRQ BASE OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p $F/ldsc/l2 $F/sumstats $F/logs/chr $F/results

one_chr () {
  set -e
  c=$1
  [[ -s "$F/ldsc/models/RMJ/RMJ.$c.annot.gz" ]] && return 0
  $RENV/bin/python -I $F/scripts/01_roadmap_annot.py "$F" $c "$PLINK" "$FRQ" > "$F/logs/chr/annot_$c.txt" 2>&1
  $LENV/bin/python $LENV/bin/ldsc.py --l2 --bfile "${PLINK}$c" --ld-wind-cm 1 \
    --annot "$F/ldsc/annot_l2/rm.$c.annot.gz" --thin-annot --print-snps "$A1/ldsc/print_snps/$c.snp" \
    --out "$F/ldsc/l2/rm.$c" > "$F/logs/chr/l2_$c.txt" 2>&1
  /usr/bin/time -v $RENV/bin/python -I $F/scripts/02_models.py "$F" $c "$H22" "$BASE" > "$F/logs/chr/models_$c.txt" 2>&1
  echo "chr$c done $(date -u +%T)"
}
export -f one_chr
seq 1 22 | xargs -P 8 -I{} bash -c 'one_chr {}'
for c in $(seq 1 22); do test -s "$F/ldsc/models/RMJ/RMJ.$c.annot.gz" || { echo "chr$c models missing" >&2; exit 1; }; done

$RENV/bin/python -I $F/scripts/03_pair_sumstats.py $A1/sumstats/ast_panukbb.sumstats.gz $H22/ldsc/sumstats/alt.sumstats.gz $F/sumstats/pair_panukbb
$RENV/bin/python -I $F/scripts/03_pair_sumstats.py $ROOT/GWAS/ldsc/munged/mvp_ast.sumstats.gz $ROOT/GWAS/ldsc/munged/mvp_alt.sumstats.gz $F/sumstats/pair_mvp
sha256sum $F/sumstats/pair_*.sumstats.gz > $F/pair_sumstats_sha256.txt

cat > $F/ldsc/config.tsv <<CFG
ldsc_python	$LENV/bin/python
ldsc_py	$LENV/bin/ldsc.py
baseline	$BASE
weights	$REF/sldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
frq	$FRQ
h22	$H22
CFG
awk -F'\t' 'NR>1 && $1 ~ /^(ast_panukbb|ast_mvp|alt_panukbb|alt_mvp|ggt_panukbb|pdff|height_panukbb)$/ {print $1"\t"$NF}' $A1/traits.tsv >> $F/ldsc/config.tsv
test $(wc -l < $F/ldsc/config.tsv) -eq 13
$RENV/bin/python -I $F/scripts/04_h2.py "$F" $F/ldsc/config.tsv 10
$RENV/bin/python -I $F/scripts/05_summarize.py "$F" "$BASE" "$H22"
$RENV/bin/python -m pip freeze > $F/pip_freeze.txt 2>/dev/null || true
echo "A1 follow-up done"
