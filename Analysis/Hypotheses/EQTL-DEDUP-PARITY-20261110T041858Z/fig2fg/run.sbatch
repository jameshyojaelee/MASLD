#!/bin/bash
#SBATCH --job-name=fig2fg_refit
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/fig2fg
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/fig2fg/run_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/fig2fg/run_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
$PY -I /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/fig2fg/merge_coloc.py
$PY -I /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/fig2fg/73_refit.py /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/fig2fg/out_loci_$SLURM_JOB_ID
$PY -m pip freeze > /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/fig2fg/out_loci_$SLURM_JOB_ID/pip_freeze.txt
