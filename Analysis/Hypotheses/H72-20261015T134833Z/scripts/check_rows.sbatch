#!/bin/bash
#SBATCH --job-name=6124
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H72-20261015T134833Z/logs/6124-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H72-20261015T134833Z/logs/6124-%j.err
# Independent re-read of two F1 rows straight from the Neale files (zcat + awk), to check the run's parser.
set -euo pipefail
H72=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H72-20261015T134833Z
S=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/MR_Data/PanUKBB/sex_stratified
O=$H72/check-20261015T1445Z
for f in $S/30730_irnt.gwas.imputed_v3.female.varorder.tsv.bgz $S/30730_irnt.gwas.imputed_v3.male.varorder.tsv.bgz \
         $H72/downloads/30730_raw.gwas.imputed_v3.female.varorder.tsv.bgz $H72/downloads/30730_raw.gwas.imputed_v3.male.varorder.tsv.bgz; do
  (zcat $f | awk -v F=$(basename $f) -F'\t' '$1=="22:24997070:A:G"{print F"\t"$1"\t"$8"\t"$9; exit}' > $O/$(basename $f).row) &
done
for f in $S/30620_irnt.gwas.imputed_v3.female.varorder.tsv.bgz $S/30620_irnt.gwas.imputed_v3.male.varorder.tsv.bgz \
         $H72/downloads/30620_raw.gwas.imputed_v3.female.varorder.tsv.bgz $H72/downloads/30620_raw.gwas.imputed_v3.male.varorder.tsv.bgz; do
  (zcat $f | awk -v F=$(basename $f) -F'\t' '$1=="1:155129562:G:A"{print F"\t"$1"\t"$8"\t"$9; exit}' > $O/$(basename $f).row) &
done
wait
cat $O/*.row > $O/rows.tsv
