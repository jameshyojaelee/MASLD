#!/bin/bash
#SBATCH --job-name=u10_rederive
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U10/rederive/logs/%j.out
# Usage: sbatch run_u10.sbatch maxT|systems
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONNOUSERSITE=1
U=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U10/rederive
E=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
MODE=$1
OUT=$U/out_${MODE}_$(date -u +%Y%m%dT%H%M%SZ)_${SLURM_JOB_ID}
mkdir -p "$OUT"
echo "mode $MODE out $OUT host $(hostname) start $(date)"
if [ "$MODE" = maxT ]; then
  cp $U/u10_maxT.R "$OUT"/
  $E/Rscript --vanilla $U/u10_maxT.R "$OUT" 2>&1 | tee "$OUT/u10_maxT.log"
else
  cp $U/u10_systems.py "$OUT"/
  $E/python3 $U/u10_systems.py "$OUT" 2>&1 | tee "$OUT/u10_systems.stdout"
  $E/python3 -m pip freeze > "$OUT/pip_freeze.txt"
fi
echo "done $(date)"
