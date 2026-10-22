#!/bin/bash
#SBATCH --job-name=9708
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9c_candidate-20261015T191932Z/build_9708_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9c_candidate-20261015T191932Z/build_9708_%j.err
export PYTHONNOUSERSITE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
D=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9c_candidate-20261015T191932Z
/usr/bin/time -v $PY $D/build_catalog_r9c.py $D/build_out > $D/build_out/build.log 2> $D/build_out/build_time_and_stderr.txt
echo build exit $?
$PY $D/check_catalog_r9c.py $D/check > $D/check/stdout.txt 2> $D/check/stderr.txt
echo check exit $?
NEGCTL=1 $PY $D/check_catalog_r9c.py $D/check/negctl > $D/check/negctl/stdout.txt 2> $D/check/negctl/stderr.txt
echo negctl exit $?
