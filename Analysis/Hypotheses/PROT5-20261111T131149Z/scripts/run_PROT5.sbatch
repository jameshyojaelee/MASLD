#!/bin/bash
#SBATCH --job-name=prot5
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PROT5-20261111T131149Z/logs/prot5_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PROT5-20261111T131149Z/logs/prot5_%j.err
# PROT5: one job, no array. Prespecification: ../PRESPEC_PROT5.md
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RUN=$ROOT/Analysis/Hypotheses/PROT5-20261111T131149Z
LOCI=$ROOT/Analysis/Hypotheses/DIG1-20261024T113951Z/out_release_vs_disease_loci/loci.tsv
OUT=$RUN/out_${SLURM_JOB_ID}
test ! -e "$OUT"; mkdir -p "$OUT"
module load PLINK/1.9
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
RS=/gpfs/commons/home/jameslee/micromamba/envs/genomicsem/bin/Rscript
echo "host $(hostname) | $(date -u) | $(plink --version | head -1)"
sha256sum -c $RUN/PRESPEC_PROT5.md.sha256 --status && echo "prespec sha256 ok" || { echo "prespec changed"; exit 1; }
sha256sum $RUN/PRESPEC_PROT5.md $RUN/scripts/* $RUN/inputs/panels.tsv $LOCI > $OUT/script_and_input_sha256.txt
cd $OUT
bash $RUN/scripts/01_download.sh $RUN
$PY -I $RUN/scripts/02_null_snps.py $ROOT $LOCI $RUN/inputs/panels.tsv $OUT
bash $RUN/scripts/03_extract.sh $RUN $OUT
$RS --vanilla $RUN/scripts/04_prot5.R $RUN $OUT
$PY -I -m pip freeze > $OUT/pip_freeze_spatial.txt
echo "done $(date -u)"
