#!/bin/bash
#SBATCH --job-name=loco_instaprism_20261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=140G
#SBATCH --time=06:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/loco_instaprism_20261118_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/loco_instaprism_20261118_%j.err

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
mkdir -p "${PROJECT_ROOT}/Analysis/Deconvolution/logs" \
         "${PROJECT_ROOT}/Analysis/Deconvolution/candidates/loco_20261118"

export R_LIBS_USER=""
export R_LIBS_SITE="${PROJECT_ROOT}/RNA-seq/.mamba/music_deconv/lib/R/library"
export LD_LIBRARY_PATH="${PROJECT_ROOT}/RNA-seq/.mamba/music_deconv/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"${PROJECT_ROOT}/RNA-seq/.mamba/music_deconv/bin/R" --vanilla --slave \
  -f "${PROJECT_ROOT}/Analysis/Deconvolution/scripts/13_run_loco_instaprism.R"
echo "LOCO_INSTAPRISM_SBATCH_DONE"
