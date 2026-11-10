#!/bin/bash
#SBATCH --job-name=eqtl_parity_E
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/logs/overlap_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/logs/overlap_%j.err
set -eo pipefail
export PYTHONNOUSERSITE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
$PY /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/code/overlap_refit.py /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/results/overlap_$SLURM_JOB_ID
$PY -m pip freeze > /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/results/overlap_$SLURM_JOB_ID/pip_freeze.txt
