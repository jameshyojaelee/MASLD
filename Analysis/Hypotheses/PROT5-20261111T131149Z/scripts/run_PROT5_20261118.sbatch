#!/bin/bash
#SBATCH --job-name=prot5_20261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/PROT5-20261111T131149Z/logs/prot5_20261118_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/PROT5-20261111T131149Z/logs/prot5_20261118_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
RUN=$ROOT/Analysis/Hypotheses/PROT5-20261111T131149Z
LOCI=$ROOT/Analysis/Hypotheses/DIG1-20261024T113951Z/out_release_vs_disease_loci_20261118/loci.tsv
OUT=$RUN/out_20261118
test ! -e "$OUT"; mkdir -p "$OUT"
module load PLINK/1.9
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
RS=/gpfs/commons/home/jameslee/micromamba/envs/genomicsem/bin/Rscript
echo "host $(hostname) | $(date -u) | $(plink --version | head -1)"
sha256sum $RUN/PRESPEC_PROT5.md $RUN/scripts/02_null_snps.py \
  $RUN/scripts/03_extract_20261118.sh $RUN/scripts/04_prot5_20261118.R \
  $RUN/scripts/run_PROT5_20261118.sbatch $RUN/inputs/panels.tsv $LOCI \
  $ROOT/Analysis/Hypotheses/ENZFAT1-20261110T204401Z/out_S2_20261118/locus_effects.tsv \
  > $OUT/script_and_input_sha256.txt
cd $OUT
$PY -I $RUN/scripts/02_null_snps.py $ROOT $LOCI $RUN/inputs/panels.tsv $OUT
bash $RUN/scripts/03_extract_20261118.sh $RUN $OUT 4
$RS --vanilla $RUN/scripts/04_prot5_20261118.R $RUN $OUT
$PY -I -m pip freeze > $OUT/pip_freeze_spatial.txt
echo "PROT5 20261118 DONE $(date -u)"
