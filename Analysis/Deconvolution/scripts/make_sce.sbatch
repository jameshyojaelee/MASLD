#!/usr/bin/env bash
#SBATCH --job-name=make_human_sce
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --partition=cpu
#SBATCH --output=RNA-seq/deconvolution/logs/make_human_sce_%j.out
#SBATCH --error=RNA-seq/deconvolution/logs/make_human_sce_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13

# Use Conda environment's libraries
export LD_LIBRARY_PATH="RNA-seq/.mamba/music_deconv/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1
export R_LIBS_USER=
export R_LIBS_SITE="RNA-seq/.mamba/music_deconv/lib/R/library"

MICROMAMBA=RNA-seq/in-house_MCD_RNAseq/micromamba
ENV_PREFIX=RNA-seq/.mamba/music_deconv

INPUT_H5AD=${INPUT_H5AD:-RNA-seq/deconvolution/reference/reference_human.consensus.h5ad}
OUTPUT_RDS=${OUTPUT_RDS:-RNA-seq/deconvolution/reference/reference_human_sce.rds}
LABEL_KEY=${LABEL_KEY:-consensus_label}

echo "Converting ${INPUT_H5AD} to ${OUTPUT_RDS} with label key ${LABEL_KEY}"
"${MICROMAMBA}" run -p "${ENV_PREFIX}" Rscript \
  RNA-seq/deconvolution/scripts/05_make_sce.R \
  "${INPUT_H5AD}" \
  "${OUTPUT_RDS}" \
  "${LABEL_KEY}"

echo "Conversion complete!"
