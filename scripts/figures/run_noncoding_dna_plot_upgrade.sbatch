#!/bin/bash
#SBATCH --job-name=ncdnaplt
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncdnaplt_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncdnaplt_%j.err

set -euo pipefail
umask 027

: "${NONCODING_DNA_INPUT:?Set NONCODING_DNA_INPUT to the immutable complete-annotation candidate}"
: "${NONCODING_DNA_PLOT_OUTPUT:?Set NONCODING_DNA_PLOT_OUTPUT to a new candidate directory}"

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"

test ! -e "${NONCODING_DNA_PLOT_OUTPUT}"

"${RNASEQ_ENV}/bin/Rscript" --vanilla scripts/figures/noncoding_dna_plot_upgrade.R \
  --input="${NONCODING_DNA_INPUT}" \
  --specificity="${PROJECT_ROOT}/GWAS/finemapping/results/gwas_atac/variant_celltype_specificity.csv" \
  --output="${NONCODING_DNA_PLOT_OUTPUT}"
