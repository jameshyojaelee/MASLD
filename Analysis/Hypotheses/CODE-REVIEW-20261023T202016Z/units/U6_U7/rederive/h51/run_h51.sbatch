#!/bin/bash
#SBATCH --job-name=U7h51
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U6_U7/rederive/h51/logs/%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
OUT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U6_U7/rederive/h51/out-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e $OUT; mkdir -p $OUT; export TMPDIR=$OUT/tmp; mkdir -p $TMPDIR
cp /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U6_U7/rederive/h51/h51_check.R $OUT/h51_check.as_run.R
echo "job $SLURM_JOB_ID host $(hostname) start $(date -u)"
/usr/bin/time -v /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U6_U7/rederive/h51/h51_check.R $OUT
echo "end $(date -u)"
