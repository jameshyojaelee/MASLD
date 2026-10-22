#!/bin/bash
#SBATCH --job-name=9705
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9b_candidate-20261015T185230Z/build_9705_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9b_candidate-20261015T185230Z/build_9705_%j.err
export PYTHONNOUSERSITE=1
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9b_candidate-20261015T185230Z
/usr/bin/time -v /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9b_candidate-20261015T185230Z/build_catalog_r9b.py /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9b_candidate-20261015T185230Z/build_out > /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9b_candidate-20261015T185230Z/build_out/build.log 2> /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9b_candidate-20261015T185230Z/build_out/build_time_and_stderr.txt
echo exit $?
