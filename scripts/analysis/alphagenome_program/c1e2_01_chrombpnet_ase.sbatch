#!/usr/bin/env bash
#SBATCH --job-name=agp_c1-endpoint2-gpu-cbp
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c1e2_cbp_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c1e2_cbp_%j.err
set -euo pipefail

# Score the 500 union allelic-imbalance sites with BOTH project-trained adult hepatocyte
# ChromBPNet families, five chromosome-group folds each. Every site is scored once, by the fold
# whose TEST group holds its chromosome, so no site is scored by a model that saw its chromosome.
# Scorer, genome, peak set, chrom sizes, batch size and random seed are the ones the deposited
# Currin calibration used, so these scores sit on the same scale as currin_scored_primary.tsv.gz.
# Nothing is written outside the package directory.

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PKG=${ROOT}/GWAS/finemapping/results/alphagenome_program/c1-endpoint2-gpu-20260914T203750Z
SCORER=${ROOT}/GWAS/finemapping/src/chrombpnet_variant_scorer/src
FASTA=/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa
BASEDIR=${ROOT}/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet

module load CUDA/11.7.0
module load cuDNN/8.4.1.50-CUDA-11.7.0
ENV=${ROOT}/.mamba/seqfunc_chrombpnet
export PATH=${ENV}/bin:${PATH}
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
export TF_FORCE_GPU_ALLOW_GROWTH=true TF_CPP_MIN_LOG_LEVEL=1

echo "host $(hostname) start $(date -u)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true
"${ENV}/bin/python" -c "import tensorflow as tf; print('tf', tf.__version__, 'gpus', tf.config.list_physical_devices('GPU'))"

for FAM in hepatocyte_5fold_v2 gse281367_hepatocyte_5fold_v2; do
  BASE=${BASEDIR}/${FAM}
  IN=${PKG}/scoring/chrombpnet/${FAM}
  OUT=${PKG}/scoring/chrombpnet/${FAM}/scores
  mkdir -p "${OUT}"
  for FOLD in 0 1 2 3 4; do
    MODEL=${BASE}/fold${FOLD}/chrombpnet_full/models/chrombpnet_nobias.h5
    test -s "${MODEL}"
    test -s "${IN}/fold${FOLD}.variants.tsv"
    if [ -s "${OUT}/fold${FOLD}.variant_scores.tsv" ]; then
      echo "[${FAM} fold${FOLD}] exists, skip"; continue
    fi
    echo "=== [${FAM} fold${FOLD}] $(date -u) n=$(wc -l < "${IN}/fold${FOLD}.variants.tsv") ==="
    "${ENV}/bin/python" "${SCORER}/variant_scoring.py" \
      -l "${IN}/fold${FOLD}.variants.tsv" \
      -g "${FASTA}" \
      -m "${MODEL}" \
      -o "${OUT}/fold${FOLD}" \
      -s "${BASE}/GRCh38.autosomes.chrom.sizes" \
      -sc chrombpnet \
      -p "${BASE}/hepatocyte.idr_pooled_summit.narrowPeak" \
      -mp 30000 -n 0 -bs 1024 -r 42 --no_hdf5
    test -s "${OUT}/fold${FOLD}.variant_scores.tsv"
    echo "[${FAM} fold${FOLD}] OK $(wc -l < "${OUT}/fold${FOLD}.variant_scores.tsv") rows"
  done
done
echo "done $(date -u)"
