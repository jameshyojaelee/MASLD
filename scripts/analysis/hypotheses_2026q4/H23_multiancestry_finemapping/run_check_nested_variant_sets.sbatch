#!/bin/bash
#SBATCH --job-name=H23_varsets
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H23-20261009T140416Z/logs/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H23-20261009T140416Z/logs/%x-%j.err
# Usage: sbatch run_check_nested_variant_sets.sbatch <fits_dir> <new_output_dir>
set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u
Rscript scripts/analysis/hypotheses_2026q4/H23_multiancestry_finemapping/check_nested_variant_sets.R "$1" "$2"
