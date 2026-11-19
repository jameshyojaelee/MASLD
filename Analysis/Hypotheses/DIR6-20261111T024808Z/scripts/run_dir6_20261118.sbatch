#!/bin/bash
#SBATCH --job-name=DIR6_20261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/DIR6-20261111T024808Z/logs/dir6_20261118_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/DIR6-20261111T024808Z/logs/dir6_20261118_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
RUN=$ROOT/Analysis/Hypotheses/DIR6-20261111T024808Z
S1=$RUN/out_23190132
OUT=$RUN/out_20261118
if [[ -e $OUT ]]; then echo "$OUT exists; not overwriting"; exit 1; fi
mkdir -p $OUT
echo "host $(hostname) | $(date -u)"
for f in rows.tsv instruments.tsv controls_variants.tsv proxies.tsv orientation_asserts.tsv pip_freeze_spatial.txt; do
  cp $S1/$f $OUT/
done
sha256sum $RUN/PRESPEC_DIR6.md $RUN/PRESPEC_DIR6_AMENDMENT1_coeffects.md \
  $RUN/scripts/02_classify_20261118.R $RUN/scripts/04_coeffects_20261118.py \
  $RUN/scripts/run_dir6_20261118.sbatch $S1/rows.tsv $S1/instruments.tsv \
  $ROOT/Analysis/Hypotheses/ENZFAT1-20261110T204401Z/out_S2_20261118/locus_effects.tsv \
  > $OUT/script_and_prespec_sha256.txt

/gpfs/commons/home/jameslee/micromamba/envs/genomicsem/bin/Rscript --vanilla \
  $RUN/scripts/02_classify_20261118.R $ROOT $OUT

for a in GCST90239658 GCST90475667 GCST90132314; do
  (cd $RUN/downloads_untrusted/$a && awk -v f=$a.h.tsv.gz '$2==f{print $1"  "$2}' harmonised_md5sum.txt | md5sum -c -)
done
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
$PY -I $RUN/scripts/04_coeffects_20261118.py $ROOT $RUN $OUT $OUT/coeffects
echo "DIR6 20261118 DONE $(date -u)"
