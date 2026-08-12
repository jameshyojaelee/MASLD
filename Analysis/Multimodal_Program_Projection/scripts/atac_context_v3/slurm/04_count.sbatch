#!/bin/bash
#SBATCH --job-name=ataccount
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_count_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_count_%j.err

set -euo pipefail
if [[ $# -ne 1 || ! $1 =~ ^GSE(244832|281367)$ ]]; then
  echo "Usage: sbatch 04_count.sbatch GSE244832|GSE281367" >&2
  exit 2
fi
COHORT=$1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PACKAGE=${ROOT}/Analysis/Multimodal_Program_Projection/scripts/atac_context_v3
SCRATCH=/scratch/${USER}/atac_context_v3_count_${COHORT}_${SLURM_JOB_ID}
MAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
mkdir -p "${SCRATCH}"
trap 'rm -rf "${SCRATCH}"' EXIT
eval "$(${MAMBA} shell hook --shell bash)"
micromamba activate snapatac2
python "${PACKAGE}/03_count_fragments.py" --cohort "${COHORT}" --tempdir "${SCRATCH}/consensus"
python "${PACKAGE}/03_count_fragments.py" --cohort "${COHORT}" --peak-space legacy_gse244832 --tempdir "${SCRATCH}/legacy"
