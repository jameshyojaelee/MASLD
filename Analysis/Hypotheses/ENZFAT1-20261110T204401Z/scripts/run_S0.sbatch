#!/bin/bash
#SBATCH --job-name=ENZFAT1_S0
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ENZFAT1-20261110T204401Z/logs/S0_%j.out
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RUN=$ROOT/Analysis/Hypotheses/ENZFAT1-20261110T204401Z
OUT=$RUN/out_S0
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
LD=/gpfs/commons/home/jameslee/micromamba/envs/ldsc/bin
REF=$ROOT/GWAS/ldsc_ref/eur_w_ld_chr
MA=$ROOT/GWAS/ldsc_ref/hm3_merge_alleles.txt
TRAITS="pdff alt_ukb alt_mvp cirr_ghouse cirr_mvp nafld_mvp nafld_finngen ast_ukb ggt_ukb"
if [[ -e $OUT ]]; then echo "$OUT exists; not overwriting"; exit 1; fi
mkdir -p $OUT/inputs $OUT/checks $OUT/munged $OUT/h2
cd $ROOT
sha256sum $RUN/scripts/00_build_input.py $RUN/scripts/01_orientation.py $RUN/scripts/run_S0.sbatch > $OUT/script_sha256.txt
M=GWAS/MR_Data
sha256sum $M/GCST90267352_PDFF_Pazoki2022.tsv.gz \
  $M/PanUKBB/biomarkers-30620-ALT-both_sexes-irnt.tsv.bgz $M/PanUKBB/biomarkers-30650-AST-both_sexes-irnt.tsv.bgz \
  $M/PanUKBB/biomarkers-30730-GGT-both_sexes-irnt.tsv.bgz \
  $M/MVP/raw/labs/MVP_R4.1000G_AGR.ALT_Mean_INT.EUR.GIA.dbGaP.txt.gz \
  $M/MVP/raw/phecodes/MVP_R4.1000G_AGR.Phe_571_51.EUR.GIA.dbGaP.txt.gz \
  $M/MVP/raw/phecodes/MVP_R4.1000G_AGR.Phe_571_5.EUR.GIA.dbGaP.txt.gz \
  $M/Ghouse_Cirrhosis/GCST90319877_harmonised_hg38.tsv.gz $M/FinnGen/finngen_R12_NAFLD.gz \
  GWAS/ldsc_ref/hm3_merge_alleles.txt GWAS/ldsc_ref/hm3_chrbp_to_rsid.tsv > $OUT/inputs_sha256.txt &
$PY -I -m pip freeze > $OUT/pip_freeze_spatial.txt
$LD/python -m pip freeze > $OUT/pip_freeze_ldsc.txt 2>/dev/null || true

# two traits at a time (2 CPUs)
echo $TRAITS | tr ' ' '\n' | xargs -P 2 -I{} sh -c "$PY -I $RUN/scripts/00_build_input.py {} $ROOT $OUT > $OUT/inputs/{}.log 2>&1"
cat $OUT/inputs/*.log
wait
$PY -I $RUN/scripts/01_orientation.py $OUT

for t in $TRAITS; do
  $LD/python $LD/munge_sumstats.py --sumstats $OUT/inputs/$t.tsv.gz --snp SNP --a1 A1 --a2 A2 \
    --signed-sumstats BETA,0 --p P --N-col N --merge-alleles $MA --chunksize 500000 --out $OUT/munged/$t \
    > $OUT/munged/$t.stdout 2>&1
  $LD/python $LD/ldsc.py --h2 $OUT/munged/$t.sumstats.gz --ref-ld-chr $REF/ --w-ld-chr $REF/ \
    --out $OUT/h2/$t > /dev/null 2>&1
done
# observed-scale h2, its z, intercept and mean chi-square per trait
printf "trait\tn_snps\th2_obs\th2_se\th2_z\tintercept\tintercept_se\tmean_chisq\tlambda_gc\n" > $OUT/h2_summary.tsv
for t in $TRAITS; do
  awk -v t=$t '
    /^After merging with regression SNP LD/ {n=$(NF-2)}
    /^Total Observed scale h2:/ {h=$5; s=$6; gsub(/[()]/,"",s)}
    /^Lambda GC:/ {l=$3} /^Mean Chi\^2:/ {m=$3}
    /^Intercept:/ {i=$2; is=$3; gsub(/[()]/,"",is)}
    END {printf "%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t%s\t%s\n", t, n, h, s, h/s, i, is, m, l}' $OUT/h2/$t.log >> $OUT/h2_summary.tsv
done
column -t $OUT/h2_summary.tsv
echo "S0 DONE"
