#!/bin/bash
#SBATCH --job-name=H23_an
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=Analysis/Hypotheses/H23-20261009T140416Z/logs/%x-%j.out
#SBATCH --error=Analysis/Hypotheses/H23-20261009T140416Z/logs/%x-%j.err
# H23 stage 2 only. Usage: sbatch run_h23_analysis_only.sbatch <fits_dir> <new_analysis_dir> [min_units]
set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u
export H23_MIN_UNITS="${3:-10}"
Rscript scripts/analysis/hypotheses_2026q4/H23_multiancestry_finemapping/h23_analyze.R "${1:?fits}" "${2:?out}"
