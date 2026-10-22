#!/bin/bash
#SBATCH --job-name=H36_sqtl_coloc
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H36-20261009T180345Z/logs/slurm-sqtl_coloc-%j.out

# H36 step: sQTL download (eQTL Catalogue r7 QTD000270, open FTP) and colocalization. One job, no array.
set -euo pipefail
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
OUT=$REPO/Analysis/Hypotheses/H36-20261009T180345Z
SC=$REPO/scripts/analysis/hypotheses_2026q4/H36_splicing
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MASLD_PROJECT_ROOT=$REPO
echo "host $(hostname)  start $(date -u +%FT%TZ)  job ${SLURM_JOB_ID:-NA}"
test -f "$OUT/PRESPEC_H36.md"
test -f "$OUT/build/junctions.tsv.gz"
test ! -e "$OUT/sqtl_coloc"
mkdir -p "$OUT/sqtl_coloc"
IN=$OUT/inputs_untrusted/QTD000270
mkdir -p "$IN"
BASE=https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000270
for f in QTD000270.cc.tsv.gz QTD000270.cc.tsv.gz.tbi; do
  if [[ ! -s "$IN/$f" ]]; then curl -fsS --retry 5 --retry-delay 30 -o "$IN/$f" "$BASE/$f"; fi
done
gzip -t "$IN/QTD000270.cc.tsv.gz"
echo "size $(stat -c %s "$IN/QTD000270.cc.tsv.gz") bytes (FTP ETag size field 0xceea788 = 217622408)"
(cd "$IN" && sha256sum QTD000270.cc.tsv.gz QTD000270.cc.tsv.gz.tbi | tee sha256.txt)
sha256sum "$OUT/PRESPEC_H36.md" "$SC"/h36_04_sqtl_coloc.R | tee "$OUT/logs/script_sha256_sqtl_coloc_${SLURM_JOB_ID:-NA}.txt"
"$RS" "$SC/h36_04_sqtl_coloc.R" "$REPO" "$OUT/build" "$IN" "$OUT/sqtl_coloc"
echo "end $(date -u +%FT%TZ)"
