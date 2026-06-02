#!/usr/bin/env bash
#SBATCH --job-name=music_summary
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=02:00:00
#SBATCH --output=RNA-seq/deconvolution/logs/%x_%j.out
#SBATCH --error=RNA-seq/deconvolution/logs/%x_%j.err

set -euo pipefail

MICROMAMBA=${MICROMAMBA:-RNA-seq/in-house_MCD_RNAseq/micromamba}
ENV_PREFIX=${ENV_PREFIX:-RNA-seq/.mamba/music_deconv}

if [[ -x "${MICROMAMBA}" ]]; then
  "${MICROMAMBA}" run -p "${ENV_PREFIX}" Rscript RNA-seq/deconvolution/scripts/09_summary_report.R \
    RNA-seq/deconvolution/results \
    RNA-seq/deconvolution/results/summary_report.md
else
  Rscript RNA-seq/deconvolution/scripts/09_summary_report.R \
    RNA-seq/deconvolution/results \
    RNA-seq/deconvolution/results/summary_report.md
fi
