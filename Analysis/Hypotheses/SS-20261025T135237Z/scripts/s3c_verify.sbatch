#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/s3c/verify_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/s3c/verify_%j.err
# SS S3c reproduction check: build S2-layout mocks from GSE174478 and GSE240729 featureCounts and run
# s3c_continuum_transport.R on each. Expected (cft-v1-20260828T141625Z): rho 0.510 n 93 p 0.117; rho 0.710 n 66 p 0.032.
set -euo pipefail
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
V=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/verify/s3c-20261025T143842Z
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
for C in GSE174478 GSE240729; do
  $RS /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts/s3c_mock_from_featurecounts.R $C $V/mock_$C
  $RS /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts/s3c_continuum_transport.R --input $V/mock_$C/mock_s2 --stage $V/mock_$C/stage.tsv \
     --out $V/out_$C --fractions $V/mock_$C/fractions.tsv --cohort $C
done
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -I -m pip freeze > $V/pip_freeze.txt 2>/dev/null || true
echo EXIT_OK
