#!/usr/bin/env bash
#SBATCH --job-name=chromBPNet
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --array=0-4%5
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/logs/bias_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/logs/bias_%A_%a.err
set -euo pipefail

module load CUDA/11.7.0
module load cuDNN/8.4.1.50-CUDA-11.7.0
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV=${ROOT}/.mamba/seqfunc_chrombpnet
export PATH=${ENV}/bin:${PATH}
BASE=${ROOT}/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2
FOLD_ID=${SLURM_ARRAY_TASK_ID}
FOLD=${BASE}/fold${FOLD_ID}
FASTA=/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa
export PYTHONNOUSERSITE=1
export TF_FORCE_GPU_ALLOW_GROWTH=true
export TMPDIR=${FOLD}/tmp_bias_full
mkdir -p "${TMPDIR}"

test -s "${FOLD}/fold${FOLD_ID}.json"
test -s "${FOLD}/hepatocyte.nonpeaks_negatives.bed"
"${ENV}/bin/python" -c 'import tensorflow as tf; assert tf.config.list_physical_devices("GPU")'
"${ENV}/bin/chrombpnet" bias pipeline \
  -ifrag "${BASE}/fragments/hepatocyte.pooled.fragments.tsv.gz" \
  -d ATAC \
  -g "${FASTA}" \
  -c "${BASE}/GRCh38.autosomes.chrom.sizes" \
  -p "${BASE}/hepatocyte.idr_pooled_summit.narrowPeak" \
  -n "${FOLD}/hepatocyte.nonpeaks_negatives.bed" \
  -fl "${FOLD}/fold${FOLD_ID}.json" \
  -o "${FOLD}/bias_full" \
  -b 0.5 \
  -e 100 \
  -es 10 \
  -s 42 \
  -bs 64 \
  --tmpdir "${TMPDIR}"
test -s "${FOLD}/bias_full/models/bias.h5"
find "${FOLD}/bias_full" -type f -printf '%P\t%s\n' | sort > "${FOLD}/bias_full/artifact_manifest.tsv"
