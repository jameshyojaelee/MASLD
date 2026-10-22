#!/bin/bash
#SBATCH --job-name=H22_l2
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=90:00:00
# H22 job B: per-chromosome annotations (03 chrom, 03 l2annot), LD scores for
# the testable annotation columns in one ldsc --l2 pass per chromosome (1 cM
# window, 1000G EUR Phase 3, --print-snps = baselineLD SNP column), then
# packing (04).
# First attempt (job 22127239) ran --l2 on all 757 columns, including 328
# untestable ones, 7 chromosomes at a time, and was killed for memory; its
# partial outputs are in ldsc/l2_all_job22127239_oom/ and are not used.
# Usage: sbatch --output=<OUT>/logs/%x_%j.out --error=<OUT>/logs/%x_%j.err run_B_ldscores.sbatch <OUT>
set -euo pipefail
OUT="${1:?out dir required}"
test -f "$OUT/ldsc/models.tsv"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
HERE=$ROOT/scripts/analysis/hypotheses_2026q4/H22_sldsc_trait_difference
RENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
LENV=/gpfs/commons/home/jameslee/micromamba/envs/ldsc
PLINK=$ROOT/data/gsmap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC.
FRQ=$ROOT/GWAS/ldsc_ref/sldsc/1000G_Phase3_frq/1000G.EUR.QC.
BASE=$ROOT/GWAS/ldsc_ref/sldsc/baselineLD_v2.2/baselineLD.
export OUT HERE RENV LENV PLINK FRQ BASE
mkdir -p "$OUT/ldsc/l2_testable" "$OUT/ldsc/print_snps" "$OUT/logs/l2"

for c in $(seq 1 22); do
  [[ -s "$OUT/ldsc/print_snps/$c.snp" ]] || zcat "${BASE}${c}.l2.ldscore.gz" | awk 'NR>1{print $2}' > "$OUT/ldsc/print_snps/$c.snp"
done

annot_one () {
  set -e
  if [[ ! -s "$OUT/ldsc/annot_stats/$1.tsv" ]]; then
    OMP_NUM_THREADS=1 $RENV/bin/python -I $HERE/03_make_annot.py chrom "$OUT" "$1" "$PLINK" "$FRQ" \
      > "$OUT/logs/l2/annot_chr$1.txt" 2>&1
  fi
  OMP_NUM_THREADS=1 $RENV/bin/python -I $HERE/03_make_annot.py l2annot "$OUT" "$1" "$PLINK" \
    > "$OUT/logs/l2/l2annot_chr$1.txt" 2>&1
}
l2_one () {
  set -e
  OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 $LENV/bin/python $LENV/bin/ldsc.py --l2 \
    --bfile "${PLINK}$1" --ld-wind-cm 1 --annot "$OUT/ldsc/annot_l2/all.$1.annot.gz" --thin-annot \
    --print-snps "$OUT/ldsc/print_snps/$1.snp" --out "$OUT/ldsc/l2_testable/all.$1" \
    > "$OUT/logs/l2/l2_chr$1.txt" 2>&1
  OMP_NUM_THREADS=1 $RENV/bin/python -I $HERE/04_pack_ldscores.py "$OUT" "$1" "$BASE" \
    > "$OUT/logs/l2/pack_chr$1.txt" 2>&1
}
export -f annot_one l2_one

seq 1 22 | xargs -P 6 -I{} bash -c 'annot_one {}'
seq 1 22 | xargs -P 6 -I{} bash -c 'l2_one {}'

for c in $(seq 1 22); do
  test -s "$OUT/ldsc/packed/l2.$c.npy" || { echo "chr$c missing packed LD scores" >&2; exit 1; }
done
echo "H22 job B done"
