#!/bin/bash
#SBATCH --job-name=code-review
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu,io
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/s3/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/s3/%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
# SS job S3d: frozen programs vs stage (PRESPEC_SS_addendum_S3d.md), with the S3e ILR beside fit. Needs S2DIR, STAGE, S3EDIR, OUT.
set -euo pipefail
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $SS/scripts/s3d_programs.R --input "$S2DIR" --stage "$STAGE" --out "$OUT" \
  --composition "$S3EDIR/ilr7.tsv"
echo EXIT_OK
