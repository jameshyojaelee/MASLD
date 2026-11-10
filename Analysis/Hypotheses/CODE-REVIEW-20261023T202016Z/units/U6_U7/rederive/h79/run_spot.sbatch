#!/bin/bash
#SBATCH --job-name=U7h79_spot
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U6_U7/rederive/h79/logs/%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
OUT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U6_U7/rederive/h79/out_spot-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e $OUT; mkdir -p $OUT; export TMPDIR=$OUT/tmp; mkdir -p $TMPDIR
cp /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U6_U7/rederive/h79/h79_resample.R $OUT/h79_resample.as_run.R
echo "job $SLURM_JOB_ID host $(hostname) start $(date -u)"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U6_U7/rederive/h79/h79_resample.R --out $OUT --cores 9 --mode spot --seed 20260923
echo "end $(date -u)"
