#!/bin/bash
#SBATCH --job-name=A1_l2
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=72:00:00
# A1 job B (adapted from H22 job B): per chromosome, annotations (03), one
# ldsc --l2 pass over every testable column (1 cM window, 1000G EUR Phase 3,
# --print-snps = baselineLD SNPs), then packing (04). Chromosomes run
# <N_PAR> at a time, largest first; chromosomes already packed are skipped.
# Usage: sbatch --output=<OUT>/logs/%x_%j.out --error=<OUT>/logs/%x_%j.err run_B_ldscores.sbatch <OUT> <N_PAR>
set -euo pipefail
OUT="${1:?out dir required}"
N_PAR="${2:?parallel chromosomes required}"
test -f "$OUT/ldsc/models.tsv"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
HERE=$OUT/scripts
RENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
LENV=/gpfs/commons/home/jameslee/micromamba/envs/ldsc
PLINK=$ROOT/data/gsmap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC.
FRQ=$ROOT/GWAS/ldsc_ref/sldsc/1000G_Phase3_frq/1000G.EUR.QC.
BASE=$ROOT/GWAS/ldsc_ref/sldsc/baselineLD_v2.2/baselineLD.
export OUT HERE RENV LENV PLINK FRQ BASE
mkdir -p "$OUT/ldsc/l2" "$OUT/logs/l2"

one_chr () {
  set -e
  c=$1
  [[ -e "$OUT/ldsc/packed/annot.$c.done" ]] && return 0
  if [[ -s "$OUT/ldsc/packed/l2.$c.npy" ]]; then
    $RENV/bin/python -I $HERE/04_pack_ldscores.py "$OUT" $c "$BASE" > "$OUT/logs/l2/pack_chr$c.txt" 2>&1
    echo "chr$c repacked $(date -u +%T)"; return 0
  fi
  /usr/bin/time -v env OMP_NUM_THREADS=1 $RENV/bin/python -I $HERE/03_make_annot.py "$OUT" $c "$PLINK" "$FRQ" \
    > "$OUT/logs/l2/annot_chr$c.txt" 2>&1
  /usr/bin/time -v env OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 $LENV/bin/python $LENV/bin/ldsc.py --l2 \
    --bfile "${PLINK}$c" --ld-wind-cm 1 --annot "$OUT/ldsc/annot_l2/all.$c.annot.gz" --thin-annot \
    --print-snps "$OUT/ldsc/print_snps/$c.snp" --out "$OUT/ldsc/l2/all.$c" > "$OUT/logs/l2/l2_chr$c.txt" 2>&1
  /usr/bin/time -v env OMP_NUM_THREADS=1 $RENV/bin/python -I $HERE/04_pack_ldscores.py "$OUT" $c "$BASE" \
    > "$OUT/logs/l2/pack_chr$c.txt" 2>&1
  echo "chr$c done $(date -u +%T)"
}
export -f one_chr

seq 1 22 | xargs -P "$N_PAR" -I{} bash -c 'one_chr {}'

for c in $(seq 1 22); do
  test -e "$OUT/ldsc/packed/annot.$c.done" || { echo "chr$c missing packed LD scores" >&2; exit 1; }
done
echo "A1 job B done"
