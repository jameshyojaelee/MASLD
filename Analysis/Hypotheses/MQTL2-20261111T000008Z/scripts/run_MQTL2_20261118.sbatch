#!/bin/bash
#SBATCH --job-name=mqtl2_20261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/MQTL2-20261111T000008Z/logs/mqtl2_20261118_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/MQTL2-20261111T000008Z/logs/mqtl2_20261118_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
RUN=$ROOT/Analysis/Hypotheses/MQTL2-20261111T000008Z
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
export PATH=$ENV:$PATH
OUT=$RUN/out_20261118
test ! -e "$OUT"; mkdir -p "$OUT"
LOCI=$RUN/inputs/alt_loci_20261118_plus_ctrl.tsv
echo "host $(hostname) | $(date -u)"
Rscript $RUN/scripts/01_make_loci_20261118.R "$LOCI"
sha256sum $RUN/scripts/01_make_loci_20261118.R $RUN/scripts/02_coloc_20261118.R \
  $RUN/scripts/03_summarise_20261118.R $RUN/scripts/run_MQTL2_20261118.sbatch \
  $LOCI $RUN/PRESPEC_MQTL2.md \
  $ROOT/Analysis/Hypotheses/DIG1-20261024T113951Z/out_release_vs_disease_loci_20261118/loci.tsv \
  $ROOT/Analysis/Hypotheses/ENZFAT1-20261110T204401Z/out_S2_20261118/locus_effects.tsv \
  > $OUT/script_and_input_sha256.txt
mkdir -p "$OUT/smoke"
MQTL2_TEST_N=1 Rscript $RUN/scripts/02_coloc_20261118.R "$OUT/smoke" "$LOCI" 2 > "$OUT/smoke/log.txt" 2>&1 || { cat "$OUT/smoke/log.txt"; exit 1; }
tail -n 10 "$OUT/smoke/log.txt"
Rscript $RUN/scripts/02_coloc_20261118.R "$OUT" "$LOCI" 16
Rscript $RUN/scripts/03_summarise_20261118.R "$OUT" "$LOCI"
echo "MQTL2 20261118 DONE $(date -u)"
