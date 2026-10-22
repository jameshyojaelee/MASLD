#!/bin/bash
#SBATCH --job-name=9030
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F2d/jobs/slurm-9030-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F2d/jobs/slurm-9030-%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 TMPDIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F2d/tmp OMP_NUM_THREADS=4
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F2d
$ENV/python scripts/04a_sites.py > out/04a_sites.log 2>&1
$ENV/Rscript scripts/05_reannotate.R > out/05_reannotate.log 2>&1
$ENV/python scripts/06_context.py > out/06_context.log 2>&1
$ENV/python scripts/04_coloc.py > out/04_coloc.log 2>&1
$ENV/pip freeze > out/pip_freeze.txt 2>/dev/null || true
echo done
