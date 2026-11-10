#!/bin/bash
#SBATCH --job-name=mprip-source-LD
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/COLOC-IDENTITY-20261107T150017Z/logs/mprip_sources_%j.out
#SBATCH --error=Analysis/Hypotheses/COLOC-IDENTITY-20261107T150017Z/logs/mprip_sources_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python Analysis/Hypotheses/COLOC-IDENTITY-20261107T150017Z/code/inspect_mprip_sources.py
