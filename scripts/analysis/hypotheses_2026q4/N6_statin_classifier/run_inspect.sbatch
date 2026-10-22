#!/bin/bash
#SBATCH --job-name=N6_inspect
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N6-20261009T135702Z/logs/inspect-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N6-20261009T135702Z/logs/inspect-%j.err
set -euo pipefail
export PATH=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin:$PATH
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
Rscript scripts/analysis/hypotheses_2026q4/N6_statin_classifier/inspect_refit_objects.R
