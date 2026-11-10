#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.err
# SS job S3b: Catalog sidecar from an S3a output (see scripts/s3b_catalog_sidecar.R header). No fitting.
#   env -u SLURM_JOB_ID -u SLURM_JOBID sbatch --export=ALL,S2DIR=<out_s2_qc>,S3ADIR=<out_s3a-...> scripts/run_s3b.sbatch
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1  # one BLAS thread per forked worker (as H51)
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
STAGE=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external/SteatoSite/fastq/PRJEB58625_fibrosis_stage_20261025T133450Z.tsv
TS=$(date -u +%Y%m%dT%H%M%SZ)
$R $SS/scripts/s3b_catalog_sidecar.R --s3a "$S3ADIR" --input "$S2DIR" --stage "$STAGE" --out "$SS/out_s3b-$TS"
