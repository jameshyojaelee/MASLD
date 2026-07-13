#!/usr/bin/env bash
#SBATCH --job-name=chromBPNet
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/chromBPNet_smoke_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/chromBPNet_smoke_%j.err
set -euo pipefail
module load CUDA/11.7.0
module load cuDNN/8.4.1.50-CUDA-11.7.0
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV=${ROOT}/.mamba/seqfunc_chrombpnet
OUT=${ROOT}/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_fold0
export PYTHONNOUSERSITE=1
export TF_FORCE_GPU_ALLOW_GROWTH=true
micromamba run -p "${ENV}" python -c 'import tensorflow as tf; g=tf.config.list_physical_devices("GPU"); print(g); assert g, "TensorFlow sees no GPU"'
# Three-epoch bias-model pilot only. Full training is gated on its QC and GPU compatibility.
micromamba run -p "${ENV}" chrombpnet bias pipeline -ifrag "${OUT}/hepatocyte.fragments.tsv.gz" -d ATAC -g /gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa -c "${OUT}/GRCh38.chrom.sizes" -p "${OUT}/hepatocyte.peaks.narrowPeak" -n "${OUT}/hepatocyte.nonpeaks_negatives.bed" -fl "${OUT}/fold0.json" -o "${OUT}/bias_smoke" -b 0.5 -e 3 -es 2 -s 42 -bs 32
test -n "$(find "${OUT}/bias_smoke" -type f -name '*.h5' -print -quit)"
