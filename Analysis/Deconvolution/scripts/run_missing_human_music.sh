#!/usr/bin/env bash
set -euo pipefail

# This script runs MuSiC deconvolution for PRJNA512027, GSE126848, and GSE167523

MICROMAMBA=${MICROMAMBA:-RNA-seq/in-house_MCD_RNAseq/micromamba}
ENV_PREFIX=${ENV_PREFIX:-RNA-seq/.mamba/music_deconv}

REF_HUMAN=${REF_HUMAN:-Analysis/Deconvolution/reference/reference_human_sce.rds}
BULK_BASE=${BULK_BASE:-Analysis/Deconvolution/bulk}
RESULTS_BASE=${RESULTS_BASE:-Analysis/Deconvolution/results}
META_BASE=${META_BASE:-Analysis/Deconvolution/metadata}

mkdir -p "${BULK_BASE}" "${RESULTS_BASE}" "${META_BASE}"

if [[ -x "${MICROMAMBA}" ]]; then
  RUN_R=("${MICROMAMBA}" run -p "${ENV_PREFIX}" Rscript)
  RUN_PY=("${MICROMAMBA}" run -p "${ENV_PREFIX}" python3)
else
  RUN_R=(Rscript)
  RUN_PY=(python3)
fi

export R_LIBS_USER=
export R_LIBS_SITE="${ENV_PREFIX}/lib/R/library"
export LD_LIBRARY_PATH="${ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

echo "=== Running Missing Human Deconvolution ==="

DATASETS=(
  "PRJNA512027|human|RNA-seq/Human/Patient_Cohorts/results/PRJNA512027/counts/featurecounts/gene_counts.txt|RNA-seq/Human/Patient_Cohorts/pipelines/custom/PRJNA512027/metadata/SraRunTable.csv"
  "GSE126848|human|RNA-seq/Human/Patient_Cohorts/results/GSE126848/counts/featurecounts/gene_counts.txt|RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE126848/metadata/SraRunTable.csv"
  "GSE167523|human|RNA-seq/Human/Patient_Cohorts/results/GSE167523/counts/featurecounts/gene_counts.txt|RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE167523/metadata/SraRunTable.csv"
)

for entry in "${DATASETS[@]}"; do
  IFS='|' read -r name species counts_path src_metadata_path <<< "${entry}"
  
  echo "--- Processing ${name} ---"
  
  metadata_path="${META_BASE}/${name}_sra_metadata.tsv"
  
  if [[ ! -f "${metadata_path}" ]]; then
    "${RUN_PY[@]}" Analysis/Deconvolution/scripts/00_convert_sra_metadata.py \
      --input "${src_metadata_path}" \
      --output "${metadata_path}"
  fi
  
  bulk_dir="${BULK_BASE}/${name}"
  result_dir="${RESULTS_BASE}/${name}"
  
  mkdir -p "${bulk_dir}" "${result_dir}"
  
  counts_tsv="${bulk_dir}/${name}_counts.tsv"
  metadata_tsv="${bulk_dir}/${name}_metadata.tsv"
  
  echo "  Running 06_prepare_bulk_counts.R..."
  "${RUN_R[@]}" Analysis/Deconvolution/scripts/06_prepare_bulk_counts.R \
    "${counts_path}" \
    "${metadata_path}" \
    "${species}" \
    "${bulk_dir}" \
    "${name}"
    
  echo "  Running 07_run_music.R..."
  "${RUN_R[@]}" Analysis/Deconvolution/scripts/07_run_music.R \
    "${REF_HUMAN}" \
    "${counts_tsv}" \
    "${metadata_tsv}" \
    "${result_dir}" \
    "${name}"
    
  echo "  Done with ${name}!"
done

echo "=== All Missing Human Datasets Deconvolved ==="
