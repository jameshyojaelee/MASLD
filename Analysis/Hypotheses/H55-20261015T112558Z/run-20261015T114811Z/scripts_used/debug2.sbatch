#!/bin/bash
#SBATCH --job-name=5502
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H55-20261015T112558Z/logs/slurm-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H55-20261015T112558Z/logs/slurm-%j.err
set -euo pipefail
DR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H55-20261015T112558Z/debug-20261015T113647Z
mkdir -p $DR/tmp $DR/scripts_used; cp /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H55-20261015T112558Z/scripts/*.R $DR/scripts_used/
export TMPDIR=$DR/tmp XDG_CACHE_HOME=$DR/tmp OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
$RS /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H55-20261015T112558Z/scripts/make_fake_morph.R $DR/fake_morph_fat.tsv
/usr/bin/time -v $RS $DR/scripts_used/h55_stats.R $DR $DR/fake_morph_fat.tsv debug 200 20 50
$RS $DR/scripts_used/h55_report.R $DR $DR/RESULT_DEBUG.md
