#!/usr/bin/env bash
#SBATCH --job-name=chromBPNet
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/gse281367_hepatocyte_5fold_v2/logs/currin_prepare_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/gse281367_hepatocyte_5fold_v2/logs/currin_prepare_%j.err
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV=${ROOT}/.mamba/seqfunc_chrombpnet
BASE=${ROOT}/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/gse281367_hepatocyte_5fold_v2
export PYTHONNOUSERSITE=1
"${ENV}/bin/python" "${ROOT}/GWAS/finemapping/src/seqfunc_v2/chrombpnet/10_prepare_currin.py" --base "${BASE}"
