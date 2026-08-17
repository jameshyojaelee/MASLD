#!/bin/bash
#SBATCH --job-name=ncdnaval
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncdnaval_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncdnaval_%j.err

set -euo pipefail
umask 027

: "${NONCODING_DNA_INPUT:?Set NONCODING_DNA_INPUT}"
: "${NONCODING_DNA_FIGURES:?Set NONCODING_DNA_FIGURES}"
: "${NONCODING_DNA_VALIDATION:?Set NONCODING_DNA_VALIDATION to a new JSON path}"

python3 scripts/manuscript/noncoding_resource/validate_noncoding_dna_candidate.py \
  --candidate "${NONCODING_DNA_INPUT}" \
  --figures "${NONCODING_DNA_FIGURES}" \
  --output "${NONCODING_DNA_VALIDATION}"

/gpfs/commons/home/jameslee/.local/bin/ruff check \
  scripts/manuscript/noncoding_resource/annotate_credible_set_vep.py \
  scripts/manuscript/noncoding_resource/build_noncoding_dna_candidate.py \
  scripts/manuscript/noncoding_resource/validate_noncoding_dna_candidate.py

/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript -e \
  'parse(file="scripts/figures/noncoding_dna_plot_upgrade.R")'

python3 -m unittest scripts.manuscript.noncoding_resource.tests.test_noncoding_genetics
python3 scripts/manuscript/validate_resource_scope.py
