#!/usr/bin/env bash
#SBATCH --job-name=music_deconv
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=2-00:00:00
#SBATCH --output=RNA-seq/deconvolution/logs/%x_%A_%a.out
#SBATCH --error=RNA-seq/deconvolution/logs/%x_%A_%a.err
#SBATCH --array=0-4

set -euo pipefail

DATASETS=(
  "inhouse_MCD"
  "GSE156918"
  "GSE205974"
  "GSE130970"
  "GSE135251"
)

DATASET_FILTER=${DATASETS[$SLURM_ARRAY_TASK_ID]}

export MICROMAMBA=${MICROMAMBA:-RNA-seq/in-house_MCD_RNAseq/micromamba}
export ENV_PREFIX=${ENV_PREFIX:-RNA-seq/.mamba/music_deconv}

export DATASET_FILTER
export RUN_QC=1
export RUN_SUMMARY=0

RNA-seq/deconvolution/scripts/run_music_deconv.sh
