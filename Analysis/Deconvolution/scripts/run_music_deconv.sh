#!/usr/bin/env bash
set -euo pipefail

MICROMAMBA=${MICROMAMBA:-RNA-seq/in-house_MCD_RNAseq/micromamba}
ENV_PREFIX=${ENV_PREFIX:-RNA-seq/.mamba/music_deconv}

REF_MOUSE=${REF_MOUSE:-RNA-seq/deconvolution/reference/reference_mouse_sce.rds}
REF_HUMAN=${REF_HUMAN:-RNA-seq/deconvolution/reference/reference_human_sce.rds}

BULK_BASE=${BULK_BASE:-RNA-seq/deconvolution/bulk}
RESULTS_BASE=${RESULTS_BASE:-RNA-seq/deconvolution/results}
META_BASE=${META_BASE:-RNA-seq/deconvolution/metadata}

DATASET_FILTER=${DATASET_FILTER:-}
RUN_QC=${RUN_QC:-1}
RUN_SUMMARY=${RUN_SUMMARY:-1}

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

mkdir -p "${BULK_BASE}" "${RESULTS_BASE}" "${META_BASE}"

if [[ -d "${ENV_PREFIX}/lib/R/library" && ! -d "${ENV_PREFIX}/lib/R/library/MuSiC" ]]; then
  echo "MuSiC not found under ${ENV_PREFIX}/lib/R/library. Run setup_env.sh or install_r_packages.R first." >&2
fi

if [[ ! -f "${META_BASE}/GSE130970_sra_metadata.tsv" ]]; then
  "${RUN_PY[@]}" RNA-seq/deconvolution/scripts/00_convert_sra_metadata.py \
    --input RNA-seq/patient_RNAseq/data/metadata/GSE130970_SraRunTable.csv \
    --output "${META_BASE}/GSE130970_sra_metadata.tsv"
fi

if [[ ! -f "${META_BASE}/GSE135251_sra_metadata.tsv" ]]; then
  "${RUN_PY[@]}" RNA-seq/deconvolution/scripts/00_convert_sra_metadata.py \
    --input RNA-seq/patient_RNAseq/data/metadata/GSE135251_SraRunTable.csv \
    --output "${META_BASE}/GSE135251_sra_metadata.tsv"
fi

DATASETS=(
  "inhouse_MCD|mouse|RNA-seq/in-house_MCD_RNAseq/counts/featurecounts/gene_counts.txt|RNA-seq/in-house_MCD_RNAseq/metadata/samples.tsv"
  "GSE156918|mouse|RNA-seq/other_MCD_RNAseq/GSE156918/counts/featurecounts/gene_counts.txt|RNA-seq/other_MCD_RNAseq/GSE156918/metadata/samples.tsv"
  "GSE205974|mouse|RNA-seq/other_MCD_RNAseq/GSE205974/counts/featurecounts/gene_counts.txt|RNA-seq/other_MCD_RNAseq/GSE205974/metadata/samples.tsv"
  "GSE130970|human|RNA-seq/patient_RNAseq/results/GSE130970/counts/gene_counts_matrix.txt|${META_BASE}/GSE130970_sra_metadata.tsv"
  "GSE135251|human|RNA-seq/patient_RNAseq/results/GSE135251/counts/gene_counts_matrix.txt|${META_BASE}/GSE135251_sra_metadata.tsv"
)

for entry in "${DATASETS[@]}"; do
  IFS='|' read -r name species counts_path metadata_path <<< "${entry}"

  if [[ -n "${DATASET_FILTER}" && "${DATASET_FILTER}" != "${name}" ]]; then
    continue
  fi

  if [[ "${species}" == "mouse" ]]; then
    ref_sce="${REF_MOUSE}"
  else
    ref_sce="${REF_HUMAN}"
  fi

  if [[ ! -f "${ref_sce}" ]]; then
    echo "Missing reference SCE: ${ref_sce}" >&2
    exit 1
  fi

  if [[ ! -f "${counts_path}" ]]; then
    echo "Missing counts file for ${name}: ${counts_path}" >&2
    continue
  fi

  if [[ ! -f "${metadata_path}" ]]; then
    echo "Missing metadata file for ${name}: ${metadata_path}" >&2
    continue
  fi

  bulk_dir="${BULK_BASE}/${name}"
  result_dir="${RESULTS_BASE}/${name}"
  qc_dir="${result_dir}/qc"

  mkdir -p "${bulk_dir}" "${result_dir}" "${qc_dir}"

  "${RUN_R[@]}" RNA-seq/deconvolution/scripts/06_prepare_bulk_counts.R \
    "${counts_path}" \
    "${metadata_path}" \
    "${species}" \
    "${bulk_dir}" \
    "${name}"

  counts_tsv="${bulk_dir}/${name}_counts.tsv"
  metadata_tsv="${bulk_dir}/${name}_metadata.tsv"

  "${RUN_R[@]}" RNA-seq/deconvolution/scripts/07_run_music.R \
    "${ref_sce}" \
    "${counts_tsv}" \
    "${metadata_tsv}" \
    "${result_dir}" \
    "${name}"

  if [[ "${RUN_QC}" == "1" ]]; then
    "${RUN_R[@]}" RNA-seq/deconvolution/scripts/08_qc_music.R \
      "${result_dir}/${name}_music_prop_weighted.tsv" \
      "${metadata_tsv}" \
      "${qc_dir}" \
      "${name}"
  fi

done

if [[ "${RUN_SUMMARY}" == "1" ]]; then
  if compgen -G "${RESULTS_BASE}/*/qc/*_celltype_summary.tsv" > /dev/null; then
    "${RUN_R[@]}" RNA-seq/deconvolution/scripts/09_summary_report.R \
      "${RESULTS_BASE}" \
      "${RESULTS_BASE}/summary_report.md"
  else
    echo "No QC summaries found under ${RESULTS_BASE}; skipping summary report." >&2
  fi
fi
