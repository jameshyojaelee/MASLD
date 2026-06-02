#!/usr/bin/env bash
#SBATCH --job-name=music_deconv
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --output=RNA-seq/deconvolution/logs/%x_%j.out
#SBATCH --error=RNA-seq/deconvolution/logs/%x_%j.err

set -euo pipefail

if [[ -z "${DATASET_FILTER:-}" ]]; then
  echo "DATASET_FILTER is required" >&2
  exit 1
fi

export MICROMAMBA=${MICROMAMBA:-RNA-seq/in-house_MCD_RNAseq/micromamba}
export ENV_PREFIX=${ENV_PREFIX:-RNA-seq/.mamba/music_deconv}
export RUN_QC=${RUN_QC:-1}
export RUN_SUMMARY=0

RNA-seq/deconvolution/scripts/run_music_deconv.sh
