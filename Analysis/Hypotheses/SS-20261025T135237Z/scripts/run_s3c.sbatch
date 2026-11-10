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
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
# SS job S3c: continuum transport on SteatoSITE (PRESPEC_SS_addendum_S3c.md). Needs S2DIR, STAGE, S3EDIR, OUT.
set -euo pipefail
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $SS/scripts/s3c_continuum_transport.R --input "$S2DIR" --stage "$STAGE" \
  --out "$OUT" --fractions "$S3EDIR/fractions.tsv" --cohort SteatoSITE
echo EXIT_OK
