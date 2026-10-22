#!/bin/bash
#SBATCH --job-name=6303
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H90-20261015T120337Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H90-20261015T120337Z/logs/%j.err
set -euo pipefail
export H90_IN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H90-20261015T120337Z/debug_synthetic/in H90_OUT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H90-20261015T120337Z/debug_synthetic/out H90_TEST=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H90-20261015T120337Z/scripts/05_fit_h90.R
