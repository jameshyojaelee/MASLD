#!/bin/bash
#SBATCH --job-name=H23_pchk
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H23-20261009T140416Z/logs/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H23-20261009T140416Z/logs/%x-%j.err
set -eo pipefail
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u
Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H23_multiancestry_finemapping/check_eqtl_pval_parse.R
