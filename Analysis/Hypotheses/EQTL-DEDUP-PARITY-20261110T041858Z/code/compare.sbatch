#!/bin/bash
#SBATCH --job-name=eqtl_parity_C
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/logs/compare_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/logs/compare_%j.err
set -eo pipefail
eval "$(micromamba shell hook -s bash)"; micromamba activate rnaseq
export RUN_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z
Rscript $RUN_DIR/code/compare.R
