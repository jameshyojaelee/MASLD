#!/bin/bash
#SBATCH --job-name=6003
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H60-20261015T112434Z/logs/slurm-%x-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 H60_DRYRUN=1 TMPDIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H60-20261015T112434Z/out-dryrun-pyonly2-20261015T120604Z/tmp MPLCONFIGDIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H60-20261015T112434Z/out-dryrun-pyonly2-20261015T120604Z/tmp
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
$ENV/python /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H60-20261015T112434Z/out-dryrun-pyonly2-20261015T120604Z/make_fake_refit.py /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H60-20261015T112434Z/out-dryrun-pyonly2-20261015T120604Z
$ENV/python /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H60-20261015T112434Z/h60_signed_twas_de.py /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H60-20261015T112434Z/out-dryrun-pyonly2-20261015T120604Z
