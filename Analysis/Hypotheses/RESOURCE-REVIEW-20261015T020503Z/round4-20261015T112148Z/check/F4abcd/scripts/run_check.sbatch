#!/bin/bash
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH -o /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F4abcd/logs/%x_%j.out
#SBATCH -e /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F4abcd/logs/%x_%j.err
CK=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F4abcd
export TMPDIR=$CK/tmp OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NCORES=8
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
for spec in $RUNS; do
  IFS=: read ARM MODE BDRAW <<< "$spec"
  export ARM MODE BDRAW OUTDIR=$CK/out/${ARM}_${MODE}_B${BDRAW}
  /usr/bin/time -v $RS $CK/scripts/check_h59.R
done
