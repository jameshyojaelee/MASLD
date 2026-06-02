#!/usr/bin/env bash
#SBATCH --job-name=celltype_annot
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/deconvolution/logs/annot_%x_%j.out
#SBATCH --error=RNA-seq/deconvolution/logs/annot_%x_%j.err

set -euo pipefail

SPECIES=${SPECIES:-mouse}
MODEL=${MODEL:-Healthy_Mouse_Liver.pkl}
MANIFEST=${MANIFEST:-Liver_Atlas/run_manifest.tsv}
OUTDIR=${OUTDIR:-RNA-seq/deconvolution/reference}
USE_GPU=${USE_GPU:-}

MICROMAMBA=${MICROMAMBA:-RNA-seq/in-house_MCD_RNAseq/micromamba}
ENV_PREFIX=${ENV_PREFIX:-RNA-seq/.mamba/music_deconv}

mkdir -p RNA-seq/deconvolution/logs

if [[ ! -x "${MICROMAMBA}" ]]; then
  echo "micromamba not found at ${MICROMAMBA}" >&2
  exit 1
fi

${MICROMAMBA} run -p "${ENV_PREFIX}" env \
  PYTHONNOUSERSITE=1 \
  LD_LIBRARY_PATH="${ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}" \
  USE_GPU="${USE_GPU}" \
  RNA-seq/deconvolution/scripts/run_annotation_pipeline.sh \
    --species "${SPECIES}" \
    --manifest "${MANIFEST}" \
    --model "${MODEL}" \
    --outdir "${OUTDIR}"
