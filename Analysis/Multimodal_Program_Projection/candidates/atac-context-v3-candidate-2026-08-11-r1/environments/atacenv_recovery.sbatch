#!/bin/bash
#SBATCH --job-name=atacenv2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_env_recovery_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_env_recovery_%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
CANDIDATE=${ROOT}/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1
PREFIX=${ROOT}/.mamba/atac-context-v3-chromvar-r46
MAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
export MAMBA_ROOT_PREFIX=${ROOT}/.mamba/atac-context-v3-micromamba-root
export CONDA_PKGS_DIRS=${MAMBA_ROOT_PREFIX}/pkgs
if [[ -e "${PREFIX}" ]]; then
  echo "Refusing to overwrite chromVAR environment: ${PREFIX}" >&2
  exit 1
fi
"${MAMBA}" create --yes --prefix "${PREFIX}" --file "${CANDIDATE}/environments/chromvar_recovery.yml"
"${MAMBA}" env export --prefix "${PREFIX}" --explicit > "${PREFIX}/chromvar.explicit.txt"
"${MAMBA}" run --prefix "${PREFIX}" Rscript -e 'stopifnot(as.character(getRversion()) >= "4.5", as.character(getRversion()) < "4.6", as.character(BiocManager::version()) == "3.22", requireNamespace("chromVAR"), requireNamespace("motifmatchr"), requireNamespace("TFBSTools"), requireNamespace("JASPAR2024"), requireNamespace("BSgenome.Hsapiens.UCSC.hg38")); sessionInfo()'
