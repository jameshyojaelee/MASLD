#!/bin/bash
#SBATCH --job-name=make_sce_nf_20261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=00:45:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/make_sce_20261118_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/make_sce_20261118_%j.err

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
mkdir -p "${PROJECT_ROOT}/Analysis/Deconvolution/logs" \
         "${PROJECT_ROOT}/Analysis/Deconvolution/candidates/loco_20261118"

PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python3 \
  "${PROJECT_ROOT}/Analysis/Deconvolution/scripts/05_extract_100k_h5ad.py"

export R_LIBS_USER=""
export R_LIBS_SITE="${PROJECT_ROOT}/RNA-seq/.mamba/music_deconv/lib/R/library"
export LD_LIBRARY_PATH="${PROJECT_ROOT}/RNA-seq/.mamba/music_deconv/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${PROJECT_ROOT}/RNA-seq/.mamba/music_deconv/bin/R" --vanilla --slave \
  -f "${PROJECT_ROOT}/Analysis/Deconvolution/scripts/05_make_sce_20261118.R"

rm -rf "${PROJECT_ROOT}/Analysis/Deconvolution/candidates/loco_20261118/tmp_sce"
echo "MAKE_SCE_20261118_DONE"
