#!/usr/bin/env bash
set -euo pipefail

# Example usage:
#   RNA-seq/deconvolution/scripts/run_annotation_pipeline.sh \
#     --species mouse \
#     --manifest Liver_Atlas/run_manifest.tsv \
#     --model <celltypist_model>

while [[ $# -gt 0 ]]; do
  case "$1" in
    --species) SPECIES="$2"; shift 2;;
    --manifest) MANIFEST="$2"; shift 2;;
    --model) MODEL="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    *) echo "Unknown arg: $1"; exit 1;;
  esac
done

MICROMAMBA=${MICROMAMBA:-RNA-seq/in-house_MCD_RNAseq/micromamba}
ENV_PREFIX=${ENV_PREFIX:-RNA-seq/.mamba/music_deconv}

if [[ -x "${MICROMAMBA}" ]]; then
  RUN_PY=("${MICROMAMBA}" run -p "${ENV_PREFIX}" python3)
else
  RUN_PY=(python3)
fi

SPECIES=${SPECIES:-mouse}
MANIFEST=${MANIFEST:-Liver_Atlas/run_manifest.tsv}
OUTDIR=${OUTDIR:-RNA-seq/deconvolution/reference}

if [[ -z "${MODEL:-}" ]]; then
  echo "CellTypist model is required (--model)." >&2
  exit 1
fi

mkdir -p "${OUTDIR}"

REF_H5AD="${OUTDIR}/reference_${SPECIES}.h5ad"
CT_H5AD="${OUTDIR}/reference_${SPECIES}.celltypist.h5ad"
SV_H5AD="${OUTDIR}/reference_${SPECIES}.scanvi.h5ad"
FINAL_H5AD="${OUTDIR}/reference_${SPECIES}.consensus.h5ad"

"${RUN_PY[@]}" RNA-seq/deconvolution/scripts/01_build_reference.py \
  --manifest "${MANIFEST}" \
  --species "${SPECIES}" \
  --output "${REF_H5AD}"

"${RUN_PY[@]}" RNA-seq/deconvolution/scripts/02_celltypist.py \
  --input "${REF_H5AD}" \
  --output "${CT_H5AD}" \
  --model "${MODEL}" \
  --majority-voting

"${RUN_PY[@]}" RNA-seq/deconvolution/scripts/03_scanvi.py \
  --input "${CT_H5AD}" \
  --output "${SV_H5AD}" \
  ${USE_GPU:+--use-gpu}

"${RUN_PY[@]}" RNA-seq/deconvolution/scripts/04_consensus_labels.py \
  --input "${SV_H5AD}" \
  --output "${FINAL_H5AD}"
