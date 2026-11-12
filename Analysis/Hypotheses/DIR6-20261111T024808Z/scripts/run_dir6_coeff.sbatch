#!/bin/bash
#SBATCH --job-name=DIR6_coeff
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIR6-20261111T024808Z/logs/dir6_coeff_%j.out
set -eo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RUN=$ROOT/Analysis/Hypotheses/DIR6-20261111T024808Z
OUT=$RUN/out_coeff_${SLURM_JOB_ID}
if [[ -e $OUT ]]; then echo "$OUT exists; not overwriting"; exit 1; fi
for a in GCST90239658 GCST90475667 GCST90132314; do
  (cd $RUN/downloads_untrusted/$a && awk -v f=$a.h.tsv.gz '$2==f{print $1"  "$2}' harmonised_md5sum.txt | md5sum -c -)
done
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
$PY -I $RUN/scripts/04_coeffects.py $ROOT $RUN $OUT
sha256sum $RUN/PRESPEC_DIR6_AMENDMENT1_coeffects.md $RUN/scripts/04_coeffects.py $RUN/scripts/run_dir6_coeff.sbatch > $OUT/script_and_prespec_sha256.txt
$PY -m pip freeze > $OUT/pip_freeze_spatial.txt 2>/dev/null || true
echo "DIR6 COEFF DONE"
