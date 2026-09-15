#!/bin/bash
#SBATCH --job-name=hac-matched-null
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --qos=nslab
#SBATCH --output=docs/technical/continuum_matched_null_smoke/null-%j.log
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export OMP_NUM_THREADS=4
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript \
  scripts/analysis/histology_anchored_continuum/translation/05_matched_null.R "$1" "$2" 4
