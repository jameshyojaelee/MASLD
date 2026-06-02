#!/usr/bin/env bash
#SBATCH --job-name=creative_plots
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/%x_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
MICROMAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
ENV_PREFIX=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/music_deconv

RESULTS_BASE="${PROJECT_ROOT}/Analysis/Deconvolution/results"
BULK_BASE="${PROJECT_ROOT}/Analysis/Deconvolution/bulk"

export R_LIBS_USER=
export R_LIBS_SITE="${ENV_PREFIX}/lib/R/library"
export LD_LIBRARY_PATH="${ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

"${MICROMAMBA}" run -p "${ENV_PREFIX}" Rscript \
  "${PROJECT_ROOT}/Analysis/Deconvolution/scripts/10_make_creative_plots.R" \
  "${RESULTS_BASE}" \
  "${BULK_BASE}" \
  "music"
