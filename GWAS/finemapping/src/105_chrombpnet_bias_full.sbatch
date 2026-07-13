#!/usr/bin/env bash
#SBATCH --job-name=chromBPNet
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=72:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/chromBPNet_bias_full_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/chromBPNet_bias_full_%j.err
set -euo pipefail

module load CUDA/11.7.0
module load cuDNN/8.4.1.50-CUDA-11.7.0

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV=${ROOT}/.mamba/seqfunc_chrombpnet
IN=${ROOT}/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_fold0
OUT=${IN}/bias_full
FASTA=/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa

export PYTHONNOUSERSITE=1
export TF_FORCE_GPU_ALLOW_GROWTH=true
export TMPDIR=${IN}/tmp_bias_full
mkdir -p "${TMPDIR}"

# The smoke-test dependency proves GPU/runtime compatibility. The production
# run uses early stopping and the package's end-to-end bias QC pipeline.
micromamba run -p "${ENV}" python -c 'import tensorflow as tf; assert tf.config.list_physical_devices("GPU")'
micromamba run -p "${ENV}" chrombpnet bias pipeline \
  -ifrag "${IN}/hepatocyte.fragments.tsv.gz" \
  -d ATAC \
  -g "${FASTA}" \
  -c "${IN}/GRCh38.chrom.sizes" \
  -p "${IN}/hepatocyte.peaks.narrowPeak" \
  -n "${IN}/hepatocyte.nonpeaks_negatives.bed" \
  -fl "${IN}/fold0.json" \
  -o "${OUT}" \
  -b 0.5 \
  -e 100 \
  -es 10 \
  -s 42 \
  -bs 64 \
  --tmpdir "${TMPDIR}"

MODEL=$(find "${OUT}" -type f -name '*.h5' -print -quit)
test -n "${MODEL}"
test -s "${MODEL}"
find "${OUT}" -type f -printf '%P\t%s\n' | sort > "${OUT}/artifact_manifest.tsv"

