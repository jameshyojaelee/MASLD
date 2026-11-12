#!/bin/bash
#SBATCH --job-name=ENZFAT1_S3r2
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ENZFAT1-20261110T204401Z/logs/S3r2_%j.out
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RUN=$ROOT/Analysis/Hypotheses/ENZFAT1-20261110T204401Z
OUT=$RUN/out_S3_r2
S1=$RUN/out_S1
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
LD=/gpfs/commons/home/jameslee/micromamba/envs/ldsc/bin
REF=$ROOT/GWAS/ldsc_ref/eur_w_ld_chr
if [[ -e $OUT ]]; then echo "$OUT exists; not overwriting"; exit 1; fi
mkdir -p $OUT/ldsc
sha256sum $RUN/scripts/30_component.py $RUN/scripts/31_usergwas_chr22.R $RUN/scripts/run_S3.sbatch $S1/ldsc_all.rds $S1/S_observed.tsv $S1/ldsc_intercepts.tsv > $OUT/script_and_input_sha256.txt
$PY -I -m pip freeze > $OUT/pip_freeze_spatial.txt
$PY -I $RUN/scripts/30_component.py $ROOT $S1 $OUT
# component heritability (gate: h2 z >= 4) and rg with its inputs
$LD/python $LD/munge_sumstats.py --sumstats $OUT/alt_not_fat_component_GRCh37.tsv.gz --snp SNP --a1 A1 --a2 A2 \
  --signed-sumstats Z,0 --p P --N-col N --merge-alleles $ROOT/GWAS/ldsc_ref/hm3_merge_alleles.txt \
  --chunksize 500000 --out $OUT/ldsc/component > $OUT/ldsc/component_munge.stdout 2>&1
$LD/python $LD/ldsc.py --h2 $OUT/ldsc/component.sumstats.gz --ref-ld-chr $REF/ --w-ld-chr $REF/ --out $OUT/ldsc/component_h2 > /dev/null 2>&1
$LD/python $LD/ldsc.py --rg $OUT/ldsc/component.sumstats.gz,$RUN/out_S0/munged/alt_ukb.sumstats.gz,$RUN/out_S0/munged/pdff.sumstats.gz \
  --ref-ld-chr $REF/ --w-ld-chr $REF/ --out $OUT/ldsc/component_rg > /dev/null 2>&1
grep -A8 "Total Observed scale h2" $OUT/ldsc/component_h2.log
grep -A4 "Summary of Genetic Correlation" $OUT/ldsc/component_rg.log
/gpfs/commons/home/jameslee/micromamba/envs/genomicsem/bin/Rscript --vanilla $RUN/scripts/31_usergwas_chr22.R $S1 $OUT 8
echo "S3 DONE"
