#!/bin/bash
#SBATCH --job-name=wgcna05b_uncorr
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/WGCNA/logs/wgcna05b_uncorr_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/WGCNA/logs/wgcna05b_uncorr_%j.err

set -eo pipefail

# T2.15 sensitivity: rerun module preservation on uncorrected CPM
# (no ComBat / removeBatchEffect with group_binary in the mod matrix).

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${MASLD_PROJECT_ROOT}/Analysis/WGCNA/scripts"
mkdir -p "${MASLD_PROJECT_ROOT}/Analysis/WGCNA/logs"

source /gpfs/commons/home/jameslee/.bashrc
micromamba activate rnaseq

cd "${MASLD_PROJECT_ROOT}"
Rscript "${SCRIPT_DIR}/05b_module_preservation_uncorrected.R"
