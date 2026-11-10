#!/bin/bash
#SBATCH --job-name=a4_ldsc_smoke
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A4_multiome_links_zonation/logs/%x_%j.log
# A4 smoke test for job D: chr22 annot + l2 + pack, AST conversion + munge, baselineLD-only h2 for ALT.
set -euo pipefail
A4=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A4_multiome_links_zonation
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
LENV=/gpfs/commons/home/jameslee/micromamba/envs/ldsc
REF=$ROOT/GWAS/ldsc_ref
PLINK=$ROOT/data/gsmap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC.
FRQ=$REF/sldsc/1000G_Phase3_frq/1000G.EUR.QC.
BASE=$REF/sldsc/baselineLD_v2.2/baselineLD.
MA=$REF/hm3_merge_alleles.txt
L=$A4/work/ldsc
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p $L/sumstats $L/print_snps $L/l2 $L/smoke $A4/logs/ldsc
( $RENV/bin/python -I $A4/scripts/05a_panukbb_to_ldsc_input.py $ROOT/GWAS/MR_Data/PanUKBB/biomarkers-30650-AST-both_sexes-irnt.tsv.bgz \
    $REF/hm3_chrbp_to_rsid.tsv $MA 399482 $L/sumstats/ast_src.tsv.gz > $A4/logs/ldsc/05a_ast.txt 2>&1 && \
  $LENV/bin/python $LENV/bin/munge_sumstats.py --sumstats $L/sumstats/ast_src.tsv.gz --snp SNP --a1 A1 --a2 A2 \
    --signed-sumstats BETA,0 --p P --N-col N --merge-alleles $MA --chunksize 500000 --out $L/sumstats/ast > /dev/null ) &
P1=$!
zcat ${BASE}22.l2.ldscore.gz | awk 'NR>1{print $2}' > $L/print_snps/22.snp
$RENV/bin/python -I $A4/scripts/05b_ldsc_annot.py annot $A4 22 $PLINK $FRQ
OMP_NUM_THREADS=2 $LENV/bin/python $LENV/bin/ldsc.py --l2 --bfile ${PLINK}22 --ld-wind-cm 1 \
  --annot $L/annot_l2/all.22.annot.gz --thin-annot --print-snps $L/print_snps/22.snp --out $L/l2/all.22 > $A4/logs/ldsc/l2_chr22.txt 2>&1
grep -E "Total time|Analysis finished" $A4/logs/ldsc/l2_chr22.txt
$RENV/bin/python -I $A4/scripts/05b_ldsc_annot.py pack $A4 22 $BASE
cat $L/annot_stats/22.tsv | head -8
$LENV/bin/python $LENV/bin/ldsc.py --h2 $ROOT/Analysis/Hypotheses/H22-20261009T015913Z/ldsc/sumstats/alt.sumstats.gz \
  --ref-ld-chr $BASE --w-ld-chr $REF/sldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --print-delete-vals \
  --overlap-annot --frqfile-chr $FRQ --print-coefficients --out $L/smoke/alt_BASE > /dev/null
grep "Total Observed scale h2" $L/smoke/alt_BASE.log
wait $P1
grep -E "SNPs remain|Mean chi|Lambda" $L/sumstats/ast.log | tail -4
echo "smoke done"
