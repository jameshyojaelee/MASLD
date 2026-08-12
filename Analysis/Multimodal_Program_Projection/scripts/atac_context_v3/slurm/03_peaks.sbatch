#!/bin/bash
#SBATCH --job-name=atacpeak
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_peaks_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_peaks_%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PACKAGE=${ROOT}/Analysis/Multimodal_Program_Projection/scripts/atac_context_v3
SCRATCH=/scratch/${USER}/atac_context_v3_peaks_${SLURM_JOB_ID}
MAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
mkdir -p "${SCRATCH}"
trap 'rm -rf "${SCRATCH}"' EXIT
eval "$(${MAMBA} shell hook --shell bash)"
micromamba activate snapatac2
python "${PACKAGE}/02_call_consensus_peaks.py" --n-jobs "${SLURM_CPUS_PER_TASK}" --tempdir "${SCRATCH}"
