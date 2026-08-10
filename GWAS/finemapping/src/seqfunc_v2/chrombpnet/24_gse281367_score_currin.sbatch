#!/usr/bin/env bash
#SBATCH --job-name=chromBPNet
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --array=0-4%5
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/gse281367_hepatocyte_5fold_v2/logs/currin_score_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/gse281367_hepatocyte_5fold_v2/logs/currin_score_%A_%a.err
set -euo pipefail

module load CUDA/11.7.0
module load cuDNN/8.4.1.50-CUDA-11.7.0
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV=${ROOT}/.mamba/seqfunc_chrombpnet
export PATH=${ENV}/bin:${PATH}
BASE=${ROOT}/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/gse281367_hepatocyte_5fold_v2
SCORER=${ROOT}/GWAS/finemapping/src/chrombpnet_variant_scorer/src
FASTA=/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa
FOLD=${SLURM_ARRAY_TASK_ID}
OUT=${BASE}/currin_calibration/scores
mkdir -p "${OUT}"
export PYTHONNOUSERSITE=1
export TF_FORCE_GPU_ALLOW_GROWTH=true
test -s "${BASE}/fold${FOLD}/chrombpnet_full/models/chrombpnet_nobias.h5"
test -s "${BASE}/currin_calibration/inputs/fold${FOLD}.variants.tsv"
"${ENV}/bin/python" "${SCORER}/variant_scoring.py" \
  -l "${BASE}/currin_calibration/inputs/fold${FOLD}.variants.tsv" -g "${FASTA}" \
  -m "${BASE}/fold${FOLD}/chrombpnet_full/models/chrombpnet_nobias.h5" \
  -o "${OUT}/fold${FOLD}" -s "${BASE}/GRCh38.autosomes.chrom.sizes" -sc chrombpnet \
  -p "${BASE}/hepatocyte.idr_pooled_summit.narrowPeak" -mp 30000 -n 0 -bs 1024 -r 42 --no_hdf5
test -s "${OUT}/fold${FOLD}.variant_scores.tsv"
