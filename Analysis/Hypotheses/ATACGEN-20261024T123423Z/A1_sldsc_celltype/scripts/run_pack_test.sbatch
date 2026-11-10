#!/bin/bash
#SBATCH --job-name=A1_packtest
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
# Repack chr22 (smoke outputs) with per-model .annot.gz, then check row counts against the frq file.
set -euo pipefail
OUT="${1:?out dir required}"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
BASE=$ROOT/GWAS/ldsc_ref/sldsc/baselineLD_v2.2/baselineLD.
FRQ=$ROOT/GWAS/ldsc_ref/sldsc/1000G_Phase3_frq/1000G.EUR.QC.22.frq
/usr/bin/time -v $RENV/bin/python -I $OUT/scripts/04_pack_ldscores.py "$OUT" 22 "$BASE"
nfrq=$(( $(wc -l < $FRQ) - 1 ))
for f in $OUT/ldsc/models/*/*.22.annot.gz; do
  n=$(( $(zcat $f | wc -l) - 1 )); [[ $n -eq $nfrq ]] || { echo "row mismatch $f $n vs $nfrq" >&2; exit 1; }
done
echo "models with chr22 annot: $(ls $OUT/ldsc/models/*/*.22.annot.gz | wc -l); rows $nfrq each"
