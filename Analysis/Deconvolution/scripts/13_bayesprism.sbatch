#!/usr/bin/env bash
#SBATCH --job-name=bayesprism_deconv
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem
#SBATCH --output=RNA-seq/deconvolution/logs/bayesprism_%j.out
#SBATCH --error=RNA-seq/deconvolution/logs/bayesprism_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13

# Environment setup
export LD_LIBRARY_PATH="RNA-seq/.mamba/music_deconv/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1
export R_LIBS_USER=
export R_LIBS_SITE="RNA-seq/.mamba/music_deconv/lib/R/library"

MICROMAMBA=RNA-seq/in-house_MCD_RNAseq/micromamba
ENV_PREFIX=RNA-seq/.mamba/music_deconv

REF_MOUSE=${REF_MOUSE:-RNA-seq/deconvolution/reference/reference_mouse_sce.rds}
REF_HUMAN=${REF_HUMAN:-RNA-seq/deconvolution/reference/reference_human_sce.rds}

BULK_BASE=RNA-seq/deconvolution/bulk
RESULTS_BASE=RNA-seq/deconvolution/results

# Datasets: name|species|counts_path
DATASETS=(
  "inhouse_MCD|mouse|${BULK_BASE}/inhouse_MCD/inhouse_MCD_counts.tsv"
  "GSE156918|mouse|${BULK_BASE}/GSE156918/GSE156918_counts.tsv"
  "GSE205974|mouse|${BULK_BASE}/GSE205974/GSE205974_counts.tsv"
  "GSE130970|human|${BULK_BASE}/GSE130970/GSE130970_counts.tsv"
  "GSE135251|human|${BULK_BASE}/GSE135251/GSE135251_counts.tsv"
)

# Run for specific dataset if provided
DATASET_FILTER=${DATASET_FILTER:-}

for entry in "${DATASETS[@]}"; do
  IFS='|' read -r name species counts_path <<< "${entry}"

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
    continue
  fi

  if [[ ! -f "${counts_path}" ]]; then
    echo "Missing counts file: ${counts_path}, skipping ${name}" >&2
    continue
  fi

  result_dir="${RESULTS_BASE}/${name}"
  mkdir -p "${result_dir}"

  echo "=== Running BayesPrism for ${name} (${species}) ==="
  "${MICROMAMBA}" run -p "${ENV_PREFIX}" Rscript \
    RNA-seq/deconvolution/scripts/12_run_bayesprism.R \
    "${ref_sce}" \
    "${counts_path}" \
    "${result_dir}" \
    "${name}" \
    "celltype"

  echo "Completed: ${name}"
done

echo "BayesPrism deconvolution pipeline complete!"
