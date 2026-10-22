#!/bin/bash
#SBATCH --job-name=5308
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H53-20261015T112655Z/jobs/5308/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H53-20261015T112655Z/jobs/5308/%j.err
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H53-20261015T112655Z
export TMPDIR=$RUN/tmp PYTHONDONTWRITEBYTECODE=1 H53_WORKERS=6
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd $RUN
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript scripts/h53_job2.R --null-debug
