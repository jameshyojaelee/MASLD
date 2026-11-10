#!/bin/bash
#SBATCH --job-name=fig2fg_render
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/fig2fg
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/fig2fg/fig_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/fig2fg/fig_%j.err
set -euo pipefail
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/fig2fg/fig2fg_refit.R /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/fig2fg/panels_$SLURM_JOB_ID
