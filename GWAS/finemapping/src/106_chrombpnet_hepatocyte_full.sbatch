#!/usr/bin/env bash
#SBATCH --job-name=chromBPNet
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/chromBPNet_hepatocyte_full_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/chromBPNet_hepatocyte_full_%j.err
set -euo pipefail

module load CUDA/11.7.0
module load cuDNN/8.4.1.50-CUDA-11.7.0

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV=${ROOT}/.mamba/seqfunc_chrombpnet
IN=${ROOT}/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_fold0
BIAS=${IN}/bias_full
OUT=${IN}/chrombpnet_full
FASTA=/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa

export PYTHONNOUSERSITE=1
export TF_FORCE_GPU_ALLOW_GROWTH=true
export TMPDIR=${IN}/tmp_chrombpnet_full
mkdir -p "${TMPDIR}"

# ChromBPNet 1.0.1's bias pipeline writes this canonical artifact. Select it
# by exact path so prediction/contribution HDF5 files can never be mistaken
# for a trained bias model.
BIAS_MODEL=${BIAS}/models/bias.h5
test -s "${BIAS_MODEL}"
printf 'bias_model\t%s\n' "${BIAS_MODEL}" > "${BIAS}/selected_model.tsv"

micromamba run -p "${ENV}" python -c 'import tensorflow as tf; assert tf.config.list_physical_devices("GPU")'
micromamba run -p "${ENV}" chrombpnet pipeline \
  -ifrag "${IN}/hepatocyte.fragments.tsv.gz" \
  -d ATAC \
  -g "${FASTA}" \
  -c "${IN}/GRCh38.chrom.sizes" \
  -p "${IN}/hepatocyte.peaks.narrowPeak" \
  -n "${IN}/hepatocyte.nonpeaks_negatives.bed" \
  -fl "${IN}/fold0.json" \
  -o "${OUT}" \
  -b "${BIAS_MODEL}" \
  -e 100 \
  -es 10 \
  -s 42 \
  -bs 64 \
  --tmpdir "${TMPDIR}"

test -s "${OUT}/models/chrombpnet.h5"
test -s "${OUT}/models/chrombpnet_nobias.h5"
find "${OUT}" -type f -printf '%P\t%s\n' | sort > "${OUT}/artifact_manifest.tsv"

# These products remain apply-only mechanistic annotations. They must not be
# consumed by canonical fine-mapping, atlas assembly, or convergence scoring.
printf 'scope\tapply_only\ncanonical_prior\tuniform_susie_unchanged\n' > "${OUT}/firewall_contract.tsv"
