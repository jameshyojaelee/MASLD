#!/bin/bash
#SBATCH --job-name=A1_smoke
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=03:00:00
# A1 smoke test (team-lead rule): steps 01 and 02 genome-wide (fast), then
# 03 annotation, ldsc --l2 and 04 packing on chr22 only. No munging.
# Usage: sbatch --output=<OUT>/logs/%x_%j.out --error=<OUT>/logs/%x_%j.err run_A0_smoke.sbatch <OUT>
set -euo pipefail
OUT="${1:?out dir required}"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
HERE=$OUT/scripts
RENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
LENV=/gpfs/commons/home/jameslee/micromamba/envs/ldsc
V3=$ROOT/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1
CORR=$ROOT/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-statistics-correction-20261022
GTF=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
CHAIN=$ROOT/data/broadaway_eqtl/hg38ToHg19.over.chain
PLINK=$ROOT/data/gsmap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC.
FRQ=$ROOT/GWAS/ldsc_ref/sldsc/1000G_Phase3_frq/1000G.EUR.QC.
BASE=$ROOT/GWAS/ldsc_ref/sldsc/baselineLD_v2.2/baselineLD.
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
mkdir -p "$OUT/logs/prep" "$OUT/ldsc/print_snps" "$OUT/ldsc/l2"
$RENV/bin/Rscript $HERE/01_peak_sets.R "$OUT" "$V3" "$CORR" "$GTF" "$CHAIN"
$RENV/bin/python -I $HERE/02_annotations_models.py "$OUT"
c=22
zcat "${BASE}${c}.l2.ldscore.gz" | awk 'NR>1{print $2}' > "$OUT/ldsc/print_snps/$c.snp"
/usr/bin/time -v $RENV/bin/python -I $HERE/03_make_annot.py "$OUT" $c "$PLINK" "$FRQ" > "$OUT/logs/prep/annot_chr$c.txt" 2>&1
/usr/bin/time -v env OMP_NUM_THREADS=2 $LENV/bin/python $LENV/bin/ldsc.py --l2 --bfile "${PLINK}$c" --ld-wind-cm 1 \
  --annot "$OUT/ldsc/annot_l2/all.$c.annot.gz" --thin-annot --print-snps "$OUT/ldsc/print_snps/$c.snp" \
  --out "$OUT/ldsc/l2/all.$c" > "$OUT/logs/prep/l2_chr$c.txt" 2>&1
/usr/bin/time -v $RENV/bin/python -I $HERE/04_pack_ldscores.py "$OUT" $c "$BASE" > "$OUT/logs/prep/pack_chr$c.txt" 2>&1
grep -h -E "Elapsed|Maximum resident" "$OUT"/logs/prep/{annot,l2,pack}_chr$c.txt
echo "A1 smoke done"
