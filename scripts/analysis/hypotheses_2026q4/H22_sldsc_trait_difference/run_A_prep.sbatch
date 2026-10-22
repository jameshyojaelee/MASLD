#!/bin/bash
#SBATCH --job-name=H22_prep
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=90:00:00
# H22 job A: donor-level peak sets (01), Pan-UKBB EUR ALT and height munging (02),
# model table (03 models).
# Usage: sbatch --output=<OUT>/logs/%x_%j.out --error=<OUT>/logs/%x_%j.err run_A_prep.sbatch <OUT>
set -euo pipefail
OUT="${1:?out dir required}"
test -f "$OUT/PRESPEC_H22.md"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
HERE=$ROOT/scripts/analysis/hypotheses_2026q4/H22_sldsc_trait_difference
RENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
LENV=/gpfs/commons/home/jameslee/micromamba/envs/ldsc
H24=$ROOT/Analysis/Hypotheses/H24-20261008T211750Z/gse296875
ENDPOINTS=$ROOT/Analysis/MASLD_Model_Benchmark/executions/gse296875-phenotype-endpoints-20260825/endpoint_lock/donor_endpoints.tsv
CHAIN=$ROOT/data/broadaway_eqtl/hg38ToHg19.over.chain
MAP=$ROOT/GWAS/ldsc_ref/hm3_chrbp_to_rsid.tsv
MA=$ROOT/GWAS/ldsc_ref/hm3_merge_alleles.txt
ALT_SRC=$ROOT/GWAS/MR_Data/PanUKBB/biomarkers-30620-ALT-both_sexes-irnt.tsv.bgz
HGT_SRC=$OUT/downloads/panukbb_height/continuous-50-both_sexes-irnt.tsv.bgz
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
mkdir -p "$OUT/ldsc/sumstats"

[[ -s "$OUT/inputs_sha256.txt" ]] || sha256sum "$H24/atac_peak_by_unit.f32" "$H24/units.tsv" "$H24/peaks.tsv" "$ENDPOINTS" "$CHAIN" \
  "$ROOT/GWAS/ldsc/munged/pdff.sumstats.gz" "$ROOT/GWAS/ldsc/munged/ghouse_cirr.sumstats.gz" \
  "$ALT_SRC" "$HGT_SRC" "$MAP" "$MA" > "$OUT/inputs_sha256.txt"

# Pan-UKBB conversion runs alongside the R step.
( $RENV/bin/python -I $HERE/02_panukbb_to_ldsc_input.py "$ALT_SRC" "$MAP" "$MA" 400822 "$OUT/ldsc/sumstats/alt_src.tsv.gz" \
    > "$OUT/logs/02_alt.txt" 2>&1 ) &
P_ALT=$!
( $RENV/bin/python -I $HERE/02_panukbb_to_ldsc_input.py "$HGT_SRC" "$MAP" "$MA" 419596 "$OUT/ldsc/sumstats/height_src.tsv.gz" \
    > "$OUT/logs/02_height.txt" 2>&1 ) &
P_HGT=$!

$RENV/bin/Rscript $HERE/01_peak_sets.R "$OUT" "$H24" "$ENDPOINTS" "$CHAIN"
wait $P_ALT
wait $P_HGT

for t in alt height; do
  $LENV/bin/python $LENV/bin/munge_sumstats.py --sumstats "$OUT/ldsc/sumstats/${t}_src.tsv.gz" \
    --snp SNP --a1 A1 --a2 A2 --signed-sumstats BETA,0 --p P --N-col N \
    --merge-alleles "$MA" --chunksize 500000 --out "$OUT/ldsc/sumstats/$t" > /dev/null
done

$RENV/bin/python -I $HERE/03_make_annot.py models "$OUT"
$RENV/bin/python -I -m pip freeze > "$OUT/pip_freeze_rnaseq.txt" 2>/dev/null || true
$LENV/bin/python -m pip freeze > "$OUT/pip_freeze_ldsc.txt" 2>/dev/null || true
echo "H22 job A done"
