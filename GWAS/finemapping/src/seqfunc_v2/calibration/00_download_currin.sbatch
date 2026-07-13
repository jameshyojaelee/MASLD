#!/usr/bin/env bash
#SBATCH --job-name=caQTL
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/caQTL_download_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/caQTL_download_%j.err
set -euo pipefail
ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
export MASLD_PROJECT_ROOT=${ROOT}
bash "${ROOT}/GWAS/finemapping/src/seqfunc_v2/calibration/00_download_currin.sh"
