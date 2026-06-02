#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --array=1-22

# SuSiE-COLOC array wrapper — one chromosome per array task
# GWAS_NAME must be set via --export in the launcher script
# Each task processes only eGenes on its chromosome (CHR_FILTER=$SLURM_ARRAY_TASK_ID)

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

module load PLINK/2.0a5.13

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

export CHR_FILTER="${SLURM_ARRAY_TASK_ID}"

echo "=== SuSiE-COLOC: ${GWAS_NAME} chr${CHR_FILTER} ==="
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start: $(date)"

cd "${BASE}"
Rscript RNA-seq/35s_susie_coloc_broadaway.R

echo "=== Done: ${GWAS_NAME} chr${CHR_FILTER} ==="
echo "End: $(date)"
