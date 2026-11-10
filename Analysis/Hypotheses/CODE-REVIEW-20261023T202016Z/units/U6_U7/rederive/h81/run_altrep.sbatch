#!/bin/bash
#SBATCH --job-name=U7h81_altrep
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=30:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U6_U7/rederive/h81/logs/%j.out
# Code review U7: H81 replicates with an alternate seed base (40260923 + r instead of 20260923 + r),
# running the producer's own replicate and summary code unchanged; outputs only in this directory.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
D=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U6_U7/rederive/h81
H81=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H81-20261015T125233Z
OUT=$D/altrep_seed40260923-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e $OUT; mkdir -p $OUT/tmp $OUT/replicates
export TMPDIR=$OUT/tmp XDG_CACHE_HOME=$OUT/tmp
# guard outputs of the producer's final run are the precondition the replicate mode checks for
cp -r $H81/run-20261015T141242Z/guard $OUT/guard
sha256sum $H81/scripts/* > $OUT/producer_script_sha256.tsv
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
echo "job $SLURM_JOB_ID host $(hostname) start $(date -u)"
$ENV/bin/Rscript --vanilla $H81/scripts/01_h81_spikein.R --mode replicates --out $OUT --cores 16 --seed 40260923 --n-rep 50
$ENV/bin/Rscript --vanilla $H81/scripts/02_h81_summary.R --out $OUT --seed 20261923
echo "end $(date -u)"
