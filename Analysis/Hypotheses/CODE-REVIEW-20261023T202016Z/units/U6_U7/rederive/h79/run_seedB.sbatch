#!/bin/bash
#SBATCH --job-name=U7h79_seedB
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=56G
#SBATCH --time=40:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U6_U7/rederive/h79/logs/%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
OUT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U6_U7/rederive/h79/out_seedB-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e $OUT; mkdir -p $OUT; export TMPDIR=$OUT/tmp; mkdir -p $TMPDIR
cp /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U6_U7/rederive/h79/h79_resample.R $OUT/h79_resample.as_run.R
echo "job $SLURM_JOB_ID host $(hostname) start $(date -u)"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U6_U7/rederive/h79/h79_resample.R --out $OUT --cores 16 --mode full --seed 20261024 --b1 200 --b2 100 --obs-c1 233 --obs-c2 570
echo "end $(date -u)"
