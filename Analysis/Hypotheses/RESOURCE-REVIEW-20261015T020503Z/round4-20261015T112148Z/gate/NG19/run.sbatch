#!/bin/bash
#SBATCH --job-name=8200
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG19/pilot_%j.log
D=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG19
export TMPDIR=$D/tmp OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4; mkdir -p $TMPDIR
/usr/bin/time -v /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $D/pilot_timing.R
