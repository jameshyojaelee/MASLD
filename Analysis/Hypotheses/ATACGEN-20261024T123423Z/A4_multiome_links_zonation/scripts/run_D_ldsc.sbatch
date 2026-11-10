#!/bin/bash
#SBATCH --job-name=a4_ldsc
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A4_multiome_links_zonation/logs/%x_%j.log
# A4 job D: Pan-UKBB EUR AST and GGT munging (05a, verbatim H22 converter), per-chromosome
# annotations, LD scores (1 cM, 1000G EUR, baselineLD SNPs), packing (05b), regressions (05c), summary (05d).
set -euo pipefail
A4=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A4_multiome_links_zonation
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
LENV=/gpfs/commons/home/jameslee/micromamba/envs/ldsc
REF=$ROOT/GWAS/ldsc_ref
PLINK=$ROOT/data/gsmap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC.
FRQ=$REF/sldsc/1000G_Phase3_frq/1000G.EUR.QC.
BASE=$REF/sldsc/baselineLD_v2.2/baselineLD.
MAP=$REF/hm3_chrbp_to_rsid.tsv
MA=$REF/hm3_merge_alleles.txt
H22=$ROOT/Analysis/Hypotheses/H22-20261009T015913Z
L=$A4/work/ldsc
export A4 RENV LENV PLINK FRQ BASE L OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$A4"
test -s $L/sets.tsv.gz
mkdir -p $L/sumstats $L/print_snps $L/l2 logs/ldsc
sha256sum $ROOT/GWAS/MR_Data/PanUKBB/biomarkers-30650-AST-both_sexes-irnt.tsv.bgz $ROOT/GWAS/MR_Data/PanUKBB/biomarkers-30730-GGT-both_sexes-irnt.tsv.bgz \
  $H22/ldsc/sumstats/alt.sumstats.gz $H22/ldsc/sumstats/height.sumstats.gz $ROOT/GWAS/ldsc/munged/pdff.sumstats.gz \
  $ROOT/GWAS/ldsc/munged/ghouse_cirr.sumstats.gz $ROOT/GWAS/ldsc/munged/finngen_nafld.sumstats.gz $MAP $MA > ldsc_input_sha256.txt

munge_one () {
  set -e
  [[ -s "$L/sumstats/$1.sumstats.gz" ]] && return 0
  $RENV/bin/python -I "$A4/scripts/05a_panukbb_to_ldsc_input.py" "$2" "$MAP" "$MA" "$3" "$L/sumstats/$1_src.tsv.gz" > "$A4/logs/ldsc/05a_$1.txt" 2>&1
  $LENV/bin/python $LENV/bin/munge_sumstats.py --sumstats "$L/sumstats/$1_src.tsv.gz" \
    --snp SNP --a1 A1 --a2 A2 --signed-sumstats BETA,0 --p P --N-col N \
    --merge-alleles "$MA" --chunksize 500000 --out "$L/sumstats/$1" > /dev/null
}
export -f munge_one; export MAP MA
( munge_one ast $ROOT/GWAS/MR_Data/PanUKBB/biomarkers-30650-AST-both_sexes-irnt.tsv.bgz 399482 ) &
P1=$!
( munge_one ggt $ROOT/GWAS/MR_Data/PanUKBB/biomarkers-30730-GGT-both_sexes-irnt.tsv.bgz 400751 ) &
P2=$!

for c in $(seq 1 22); do
  [[ -s "$L/print_snps/$c.snp" ]] || zcat "${BASE}${c}.l2.ldscore.gz" | awk 'NR>1{print $2}' > "$L/print_snps/$c.snp"
done
annot_one () {
  set -e
  $RENV/bin/python -I $A4/scripts/05b_ldsc_annot.py annot "$A4" "$1" "$PLINK" "$FRQ" > "$A4/logs/ldsc/annot_chr$1.txt" 2>&1
}
l2_one () {
  set -e
  if [[ ! -s "$L/l2/all.$1.l2.ldscore.gz" ]]; then
    OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 $LENV/bin/python $LENV/bin/ldsc.py --l2 \
      --bfile "${PLINK}$1" --ld-wind-cm 1 --annot "$L/annot_l2/all.$1.annot.gz" --thin-annot \
      --print-snps "$L/print_snps/$1.snp" --out "$L/l2/all.$1" > "$A4/logs/ldsc/l2_chr$1.txt" 2>&1
  fi
  $RENV/bin/python -I $A4/scripts/05b_ldsc_annot.py pack "$A4" "$1" "$BASE" > "$A4/logs/ldsc/pack_chr$1.txt" 2>&1
}
export -f annot_one l2_one
seq 1 22 | xargs -P 8 -I{} bash -c 'annot_one {}'
seq 1 22 | xargs -P 6 -I{} bash -c 'l2_one {}'
for c in $(seq 1 22); do test -s "$L/packed/l2.$c.npy" || { echo "chr$c missing packed LD scores" >&2; exit 1; }; done
wait $P1; wait $P2
cat > $L/config.tsv <<CFG
ldsc_python	$LENV/bin/python
ldsc_py	$LENV/bin/ldsc.py
baseline	$BASE
weights	$REF/sldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
frq	$FRQ
alt	$H22/ldsc/sumstats/alt.sumstats.gz
ast	$L/sumstats/ast.sumstats.gz
ggt	$L/sumstats/ggt.sumstats.gz
pdff	$ROOT/GWAS/ldsc/munged/pdff.sumstats.gz
cirr	$ROOT/GWAS/ldsc/munged/ghouse_cirr.sumstats.gz
nafld	$ROOT/GWAS/ldsc/munged/finngen_nafld.sumstats.gz
height	$H22/ldsc/sumstats/height.sumstats.gz
CFG
$RENV/bin/python -I $A4/scripts/05c_ldsc_h2.py "$A4" $L/config.tsv 15
$RENV/bin/python -I $A4/scripts/05d_ldsc_summarize.py "$A4" "$BASE"
$LENV/bin/python -m pip freeze > pip_freeze_ldsc.txt 2>/dev/null || true
echo "A4 job D done"
