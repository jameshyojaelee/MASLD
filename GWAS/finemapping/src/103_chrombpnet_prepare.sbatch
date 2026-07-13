#!/usr/bin/env bash
#SBATCH --job-name=chromBPNet
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/chromBPNet_prep_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/chromBPNet_prep_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV=${ROOT}/.mamba/seqfunc_chrombpnet
OUT=${ROOT}/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_fold0
export PYTHONNOUSERSITE=1
python3 "${ROOT}/GWAS/finemapping/src/103_chrombpnet_hepatocyte_prepare.py"
CBP="micromamba run -p ${ENV} chrombpnet"
${CBP} prep splits -op "${OUT}/fold0" -c "${OUT}/GRCh38.chrom.sizes" -tcr chr1 chr8 -vcr chr2 chr9
${CBP} prep nonpeaks -g /gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa -o "${OUT}/hepatocyte.nonpeaks" -p "${OUT}/hepatocyte.peaks.narrowPeak" -c "${OUT}/GRCh38.chrom.sizes" -fl "${OUT}/fold0.json" -s 42
test -s "${OUT}/hepatocyte.nonpeaks_negatives.bed"
