#!/bin/bash
#SBATCH --job-name=A1_prep
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=72:00:00
# A1 job A: peak sets (01), annotation catalogue + seeded nulls + model table (02),
# munging of Pan-UKBB AST and GGT (EUR columns), Yengo 2018 BMI and GLGC TG
# (trans-ancestry, descriptive), then a timing test of annotation + LD scores on chr22.
# Usage: sbatch --output=<OUT>/logs/%x_%j.out --error=<OUT>/logs/%x_%j.err run_A_prep.sbatch <OUT>
set -euo pipefail
OUT="${1:?out dir required}"
test -f "$OUT/PRESPEC.md"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
HERE=$OUT/scripts
RENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
LENV=/gpfs/commons/home/jameslee/micromamba/envs/ldsc
V3=$ROOT/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1
CORR=$ROOT/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-statistics-correction-20261022
GTF=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
CHAIN=$ROOT/data/broadaway_eqtl/hg38ToHg19.over.chain
MAP=$ROOT/GWAS/ldsc_ref/hm3_chrbp_to_rsid.tsv
MA=$ROOT/GWAS/ldsc_ref/hm3_merge_alleles.txt
PUKB=$ROOT/GWAS/MR_Data/PanUKBB
BMI_SRC=$ROOT/GWAS/MR_Data/Yengo_2018_BMI_GWAS.txt.gz
TG_SRC=$ROOT/GWAS/MR_Data/GLGC/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_logTG_INV_ALL_with_N_1.gz
PLINK=$ROOT/data/gsmap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC.
FRQ=$ROOT/GWAS/ldsc_ref/sldsc/1000G_Phase3_frq/1000G.EUR.QC.
BASE=$ROOT/GWAS/ldsc_ref/sldsc/baselineLD_v2.2/baselineLD.
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
SS=$OUT/sumstats
mkdir -p "$SS" "$OUT/logs/prep" "$OUT/ldsc/print_snps" "$OUT/ldsc/l2"

[[ -s "$OUT/inputs_sha256.txt" ]] || sha256sum "$V3/consensus_peak_manifest.tsv" "$CORR/da/da_peak_results.tsv.gz" \
  "$V3/da/da_model_qc.tsv" "$V3"/counts/GSE244832/{hepatocyte,stellate}.{mtx.gz,donors.tsv,peaks.tsv} \
  "$V3"/counts/GSE281367/{hepatocyte,stellate}.{mtx.gz,donors.tsv,peaks.tsv} "$GTF" "$CHAIN" "$MAP" "$MA" \
  "$PUKB/biomarkers-30650-AST-both_sexes-irnt.tsv.bgz" "$PUKB/biomarkers-30730-GGT-both_sexes-irnt.tsv.bgz" \
  "$BMI_SRC" "$TG_SRC" > "$OUT/inputs_sha256.txt"

munge () {  # <name> <src> <extra munge args...>
  local name=$1 src=$2; shift 2
  [[ -s "$SS/$name.sumstats.gz" ]] && return 0
  $LENV/bin/python $LENV/bin/munge_sumstats.py --sumstats "$src" --merge-alleles "$MA" \
    --chunksize 500000 --out "$SS/$name" "$@" > /dev/null
}
(
  for t in "30650-AST 399482 ast_panukbb" "30730-GGT 400751 ggt_panukbb"; do
    set -- $t
    [[ -s "$SS/${3}_src.tsv.gz" ]] || $RENV/bin/python -I $HERE/00_panukbb_to_ldsc_input.py \
      "$PUKB/biomarkers-$1-both_sexes-irnt.tsv.bgz" "$MAP" "$MA" "$2" "$SS/${3}_src.tsv.gz" > "$OUT/logs/prep/00_$3.txt" 2>&1
    munge "$3" "$SS/${3}_src.tsv.gz" --snp SNP --a1 A1 --a2 A2 --signed-sumstats BETA,0 --p P --N-col N
  done
) &
P1=$!
(
  munge bmi_yengo "$BMI_SRC" --snp SNP --a1 Tested_Allele --a2 Other_Allele --signed-sumstats BETA,0 --p P --N-col N
  # 1,668 GLGC P values are below 1e-308 and fail munge's numeric check; drop P < 1e-300
  # (chi2 > ~1370, at or above ldsc's default chisq-max of max(80, 0.001 N)).
  if [[ ! -s "$SS/tg_glgc_multi_src_p1e300.tsv.gz" ]]; then
    zcat "$TG_SRC" | awk -F'\t' 'BEGIN{OFS="\t"; print "SNP","A1","A2","BETA","P","N"}
      NR>1 && $1 ~ /^rs/ && $16+0 >= 1e-300 {print $1,$5,$4,$14,$16,$6}' | gzip -1 > "$SS/tg_glgc_multi_src_p1e300.tsv.gz"
  fi
  munge tg_glgc_multi "$SS/tg_glgc_multi_src_p1e300.tsv.gz" --snp SNP --a1 A1 --a2 A2 --signed-sumstats BETA,0 --p P --N-col N
) &
P2=$!

[[ -s "$OUT/peaks/peaks.tsv" ]] || $RENV/bin/Rscript $HERE/01_peak_sets.R "$OUT" "$V3" "$CORR" "$GTF" "$CHAIN"
[[ -s "$OUT/ldsc/models.tsv" ]] || $RENV/bin/python -I $HERE/02_annotations_models.py "$OUT"

for c in $(seq 1 22); do
  [[ -s "$OUT/ldsc/print_snps/$c.snp" ]] || zcat "${BASE}${c}.l2.ldscore.gz" | awk 'NR>1{print $2}' > "$OUT/ldsc/print_snps/$c.snp"
done

# timing test: chromosome 22, all testable columns (skipped when the smoke test already packed chr22)
c=22
if [[ ! -s "$OUT/ldsc/packed/l2.$c.npy" ]]; then
/usr/bin/time -v $RENV/bin/python -I $HERE/03_make_annot.py "$OUT" $c "$PLINK" "$FRQ" > "$OUT/logs/prep/annot_chr$c.txt" 2>&1
/usr/bin/time -v env OMP_NUM_THREADS=2 $LENV/bin/python $LENV/bin/ldsc.py --l2 --bfile "${PLINK}$c" --ld-wind-cm 1 \
  --annot "$OUT/ldsc/annot_l2/all.$c.annot.gz" --thin-annot --print-snps "$OUT/ldsc/print_snps/$c.snp" \
  --out "$OUT/ldsc/l2/all.$c" > "$OUT/logs/prep/l2_chr$c.txt" 2>&1
/usr/bin/time -v $RENV/bin/python -I $HERE/04_pack_ldscores.py "$OUT" $c "$BASE" > "$OUT/logs/prep/pack_chr$c.txt" 2>&1
grep -h -E "Elapsed|Maximum resident" "$OUT"/logs/prep/{annot,l2,pack}_chr$c.txt
fi

wait $P1
wait $P2
for f in "$SS"/*.log; do echo "$f: $(grep -E 'Writing summary statistics|Mean chi' "$f" | tr '\n' ' ')"; done
$RENV/bin/python -I -m pip freeze > "$OUT/pip_freeze_rnaseq.txt" 2>/dev/null || true
$LENV/bin/python -m pip freeze > "$OUT/pip_freeze_ldsc.txt" 2>/dev/null || true
echo "A1 job A done"
