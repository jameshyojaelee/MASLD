#!/usr/bin/env bash
#SBATCH --job-name=liver_music_missing
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --time=12:00:00
#SBATCH --output=Analysis/Deconvolution/logs/music_missing_%a_%j.out
#SBATCH --error=Analysis/Deconvolution/logs/music_missing_%a_%j.err
#SBATCH --array=0-1

set -euo pipefail

# This script runs MuSiC deconvolution for GSE126848 and GSE167523

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

MICROMAMBA=${MICROMAMBA:-RNA-seq/in-house_MCD_RNAseq/micromamba}
ENV_PREFIX=${ENV_PREFIX:-RNA-seq/.mamba/music_deconv}

REF_HUMAN=${REF_HUMAN:-Analysis/Deconvolution/reference/reference_human_sce.rds}
BULK_BASE=${BULK_BASE:-Analysis/Deconvolution/bulk}
RESULTS_BASE=${RESULTS_BASE:-Analysis/Deconvolution/results}
META_BASE=${META_BASE:-Analysis/Deconvolution/metadata}

mkdir -p "${BULK_BASE}" "${RESULTS_BASE}" "${META_BASE}"

export R_LIBS_USER=
export R_LIBS_SITE="${ENV_PREFIX}/lib/R/library"
export LD_LIBRARY_PATH="${ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

DATASETS=(
  "GSE126848|human|RNA-seq/Human/Patient_Cohorts/results/GSE126848/counts/featurecounts/gene_counts.txt|RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE126848/metadata/SraRunTable.csv"
  "GSE167523|human|RNA-seq/Human/Patient_Cohorts/results/GSE167523/counts/featurecounts/gene_counts.txt|RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE167523/metadata/SraRunTable.csv"
)

entry="${DATASETS[$SLURM_ARRAY_TASK_ID]}"

IFS='|' read -r name species counts_path src_metadata_path <<< "${entry}"

echo "=== Running Missing Human Deconvolution ==="
echo "--- Processing ${name} ---"

metadata_path="${META_BASE}/${name}_sra_metadata.tsv"

if [[ ! -f "${metadata_path}" ]]; then
  python3 Analysis/Deconvolution/scripts/00_convert_sra_metadata.py \
    --input "${src_metadata_path}" \
    --output "${metadata_path}"
fi

bulk_dir="${BULK_BASE}/${name}"
result_dir="${RESULTS_BASE}/${name}"

mkdir -p "${bulk_dir}" "${result_dir}"

counts_tsv="${bulk_dir}/${name}_counts.tsv"
metadata_tsv="${bulk_dir}/${name}_metadata.tsv"

echo "  Running 06_prepare_bulk_counts.R..."
Rscript Analysis/Deconvolution/scripts/06_prepare_bulk_counts.R \
  "${counts_path}" \
  "${metadata_path}" \
  "${species}" \
  "${bulk_dir}" \
  "${name}"
  
echo "  Running 07_run_music.R..."
Rscript Analysis/Deconvolution/scripts/07_run_music.R \
  "${REF_HUMAN}" \
  "${counts_tsv}" \
  "${metadata_tsv}" \
  "${result_dir}" \
  "${name}"
  
echo "  Done with ${name}!"
