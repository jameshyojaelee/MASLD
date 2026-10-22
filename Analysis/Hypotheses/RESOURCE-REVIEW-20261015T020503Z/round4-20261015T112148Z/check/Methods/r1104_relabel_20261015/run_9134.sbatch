#!/bin/bash
#SBATCH --job-name=9134
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/Methods/r1104_relabel_20261015/slurm_9134_%j.log
export TMPDIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/Methods/r1104_relabel_20261015/tmp3; mkdir -p $TMPDIR
R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
LC_ALL=C.UTF-8 TAG=Cutf8 $R /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/Methods/r1104_relabel_20261015/identity_tables.R /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/Methods/r1104_relabel_20261015/out_identity_tables
LC_ALL=en_US.UTF-8 TAG=enUS $R /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/Methods/r1104_relabel_20261015/identity_tables.R /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/Methods/r1104_relabel_20261015/out_identity_tables
