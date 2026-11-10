#!/bin/bash
#SBATCH --job-name=code-review
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu,io
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/s3/helper_%j.out
# SS S3a parallel block helper (D-S3a-1). Submit with --export=ALL,BLOCKS=<tag:block,...>
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $SS/scripts/helper_s3a_blocks.R \
  $SS/out_s3a-20261026T234211Z $SS/out_s3a_helpers/job_${SLURM_JOB_ID} "$BLOCKS"
echo EXIT_OK
