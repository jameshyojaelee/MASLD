#!/bin/bash
#SBATCH --job-name=hac-translation
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --qos=nslab
#SBATCH --output=docs/technical/continuum_translation_20260906T135538Z/build-%j.log
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export OMP_NUM_THREADS=4
R_BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
"$R_BIN" scripts/analysis/histology_anchored_continuum/translation/test_translation.R
"$R_BIN" scripts/analysis/histology_anchored_continuum/translation/01_build.R "$1"
"$R_BIN" scripts/analysis/histology_anchored_continuum/translation/02_paired_protein.R "$1"
"$R_BIN" scripts/analysis/histology_anchored_continuum/translation/03_integrate_render.R "$1" integrated
"$R_BIN" scripts/analysis/histology_anchored_continuum/translation/04_verify.R "$1" integrated
python3 scripts/manuscript/validate_resource_scope.py
