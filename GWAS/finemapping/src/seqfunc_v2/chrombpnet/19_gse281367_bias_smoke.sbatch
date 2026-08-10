#!/usr/bin/env bash
#SBATCH --job-name=chromBPNet
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/gse281367_hepatocyte_5fold_v2/logs/bias_smoke_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/gse281367_hepatocyte_5fold_v2/logs/bias_smoke_%j.err
set -euo pipefail

module load CUDA/11.7.0
module load cuDNN/8.4.1.50-CUDA-11.7.0
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV=${ROOT}/.mamba/seqfunc_chrombpnet
export PATH=${ENV}/bin:${PATH}
BASE=${ROOT}/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/gse281367_hepatocyte_5fold_v2
FOLD=${BASE}/fold0
FASTA=/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa
export PYTHONNOUSERSITE=1
export TF_FORCE_GPU_ALLOW_GROWTH=true
export TMPDIR=${FOLD}/tmp_bias_smoke
mkdir -p "${TMPDIR}"
test -s "${BASE}/fragments/hepatocyte.pooled.fragments.tsv.gz"
test -s "${BASE}/hepatocyte.idr_pooled_summit.narrowPeak"
test -s "${FOLD}/fold0.json"
test -s "${FOLD}/hepatocyte.nonpeaks_negatives.bed"
"${ENV}/bin/python" -c 'import tensorflow as tf; g=tf.config.list_physical_devices("GPU"); print(g); assert g'
"${ENV}/bin/chrombpnet" bias pipeline \
  -ifrag "${BASE}/fragments/hepatocyte.pooled.fragments.tsv.gz" -d ATAC -g "${FASTA}" \
  -c "${BASE}/GRCh38.autosomes.chrom.sizes" -p "${BASE}/hepatocyte.idr_pooled_summit.narrowPeak" \
  -n "${FOLD}/hepatocyte.nonpeaks_negatives.bed" -fl "${FOLD}/fold0.json" \
  -o "${FOLD}/bias_smoke" -b 0.5 -e 3 -es 2 -s 42 -bs 64 --tmpdir "${TMPDIR}"
MODEL=$(find "${FOLD}/bias_smoke" -type f -name '*.h5' -print -quit)
test -n "${MODEL}"
test -s "${MODEL}"
MANIFEST_TMP=${FOLD}/bias_smoke.artifact_manifest.tsv.tmp
find "${FOLD}/bias_smoke" -type f ! -name 'artifact_manifest.tsv' -printf '%P\t%s\n' \
  | sort > "${MANIFEST_TMP}"
mv "${MANIFEST_TMP}" "${FOLD}/bias_smoke/artifact_manifest.tsv"
