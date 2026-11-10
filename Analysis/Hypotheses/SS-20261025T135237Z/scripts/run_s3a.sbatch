#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.err
# SS job S3a: F-SS1, F-SS2, F-SS3 on the S2 output (see scripts/s3a_replication.R header). 1,000 permutations x 3 fits
# plus 400 bootstrap fits. Requires PRESPEC_SS.sha256. Submit with:
#   env -u SLURM_JOB_ID -u SLURM_JOBID sbatch --export=ALL,S2DIR=<out_s2_qc>,ILR=<S3e ilr tsv or none> scripts/run_s3a.sbatch
# Restart after a timeout (reuses finished checkpoint blocks; same numbers as an uninterrupted run):
#   ... sbatch --export=ALL,S2DIR=<same>,ILR=<same>,S3A_OUT=<existing out_s3a-...> scripts/run_s3a.sbatch
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1  # one BLAS thread per forked worker (as H51)
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
STAGE=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external/SteatoSite/fastq/PRJEB58625_fibrosis_stage_20261025T133450Z.tsv
if [ -n "${S3A_OUT:-}" ]; then OUT=$S3A_OUT; RESUME=--resume; else OUT=$SS/out_s3a-$(date -u +%Y%m%dT%H%M%SZ); RESUME=; fi
$R $SS/scripts/s3a_replication.R --input "$S2DIR" --stage "$STAGE" --out "$OUT" --ilr "${ILR:-none}" --block-size 48 $RESUME
echo "S3A_OUT=$OUT"
