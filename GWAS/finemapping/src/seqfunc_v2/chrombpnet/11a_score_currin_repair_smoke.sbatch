#!/usr/bin/env bash
#SBATCH --job-name=chromBPNet
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/logs/currin_repair_smoke_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/logs/currin_repair_smoke_%j.err
set -euo pipefail

module load CUDA/11.7.0
module load cuDNN/8.4.1.50-CUDA-11.7.0
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV=${ROOT}/.mamba/seqfunc_chrombpnet
export PATH=${ENV}/bin:${PATH}
BASE=${ROOT}/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2
SCORER=${ROOT}/GWAS/finemapping/src/chrombpnet_variant_scorer/src
FASTA=/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa
OUT=${BASE}/currin_calibration/repair_smoke_20260714
mkdir -p "${OUT}"
sed -n '1p' "${BASE}/currin_calibration/inputs/fold0.variants.tsv" > "${OUT}/fold0.one_variant.tsv"
export PYTHONNOUSERSITE=1
export TF_FORCE_GPU_ALLOW_GROWTH=true

"${ENV}/bin/python" "${SCORER}/variant_scoring.py" \
  -l "${OUT}/fold0.one_variant.tsv" \
  -g "${FASTA}" \
  -m "${BASE}/fold0/chrombpnet_full/models/chrombpnet_nobias.h5" \
  -o "${OUT}/fold0" \
  -s "${BASE}/GRCh38.autosomes.chrom.sizes" \
  -sc chrombpnet \
  -p "${BASE}/hepatocyte.idr_pooled_summit.narrowPeak" \
  -mp 100 \
  -n 0 \
  -bs 1 \
  -r 42 \
  --no_hdf5
test -s "${OUT}/fold0.variant_scores.tsv"
