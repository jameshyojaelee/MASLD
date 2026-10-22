#!/bin/bash
#SBATCH --job-name=5500
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H55-20261015T112558Z/slurm-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H55-20261015T112558Z/slurm-%j.err
set -euo pipefail
D=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H55-20261015T112558Z
RUN=$D/run-${RUNSTAMP}
P=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3111
SL=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external/gtex_liver_slides/series
mkdir $RUN; mkdir -p $RUN/morph $RUN/tmp
mkdir $RUN/scripts_used; cp $D/scripts/*.R $D/scripts/*.sbatch $D/scripts/*.sha256 $RUN/scripts_used/
export TMPDIR=$RUN/tmp XDG_CACHE_HOME=$RUN/tmp MPLCONFIGDIR=$RUN/tmp PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
echo "[$(date -u +%FT%TZ)] step 1: input hashes"
sha256sum -c $D/scripts/prespec.sha256 || { echo "PRESPEC HASH MISMATCH; refusing"; exit 2; }
sha256sum -c $D/scripts/inputs.sha256 || { echo "INPUT HASH MISMATCH; refusing"; exit 2; }
NF=$(find $SL -type f | wc -l); NB=$(find $SL -type f -printf '%s\n' | awk '{s+=$1} END{printf "%d", s}')
NEWEST=$(find $SL -type f -printf '%TY-%Tm-%TdT%TH:%TM:%TS\n' | sort | tail -1 | cut -c1-19); NS=$(ls $SL | wc -l)
echo "slides: series $NS files $NF bytes $NB newest $NEWEST"
[ "$NS" = 224 ] && [ "$NF" = 1127 ] && [ "$NB" = 101783399016 ] && [ "$NEWEST" = 2026-09-23T20:49:50 ] || { echo "SLIDE INVENTORY MISMATCH; refusing"; exit 2; }
(cd $SL && find . -type f | sort | xargs -P 4 -n 20 sha256sum) | sort -k2 > $RUN/slides_sha256.txt
echo "[$(date -u +%FT%TZ)] slide files hashed: $(wc -l < $RUN/slides_sha256.txt)"
echo "INPUT_HASH_CHECKS_PASSED"
PY=/gpfs/commons/home/jameslee/micromamba/envs/he_molecular/bin/python
$PY -m pip freeze > $RUN/pip_freeze_he_molecular.txt 2>/dev/null || true
if diff -q $RUN/pip_freeze_he_molecular.txt $P/pip_freeze.txt >/dev/null; then echo "he_molecular pip freeze identical to pilot"; else echo "he_molecular pip freeze DIFFERS from pilot (see diff)"; diff $RUN/pip_freeze_he_molecular.txt $P/pip_freeze.txt > $RUN/pip_freeze_diff_vs_pilot.txt || true; fi
echo "[$(date -u +%FT%TZ)] step 2: morphometry"
cp $P/morph.py $P/order.tsv $RUN/morph/
cd $RUN/morph
echo "3c147f0c17ebff0a6c273e7cf1fd9e22f863ecf1502e63a4b68207a864324d14  morph.py" | sha256sum -c
echo "dc5a02c935b700470ea8e30573ae6f480e4965779c8c2bfeb329f6ed1353a8e2  order.tsv" | sha256sum -c
echo "7b249f92b59280434a6deeb965ff877ae9742a4ec1f29b9d21831fb2ded7ef4a  $P/morph_fat.tsv" | sha256sum -c
$PY morph.py 1000000
echo "[$(date -u +%FT%TZ)] step 3: statistics"
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
/usr/bin/time -v $RS $RUN/scripts_used/h55_stats.R $RUN $RUN/morph/morph_fat.tsv final 10000 2000 1000
$RS $RUN/scripts_used/h55_report.R $RUN $D/RESULT_H55.md
echo "[$(date -u +%FT%TZ)] done"
